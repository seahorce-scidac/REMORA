#ifndef REMORA_DateClock_H
#define REMORA_DateClock_H

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

/** \file REMORA_DateClock.H
 *
 * Model time to calendar date, ported from ROMS `dateclock.F`.
 *
 * The model clock is elapsed time since a reference date, and `remora.time_ref`
 * both sets that date and picks the calendar, exactly as ROMS `TIME_REF` does:
 *
 * | `time_ref`     | calendar              | epoch                 | year length |
 * |----------------|-----------------------|-----------------------|-------------|
 * | `yyyymmdd.dd`  | proleptic Gregorian   | the date given        | 365.2425 d  |
 * | `0`            | proleptic Gregorian   | 0001-01-01 00:00:00   | 365.2425 d  |
 * | `-1`           | 360_day               | 0000-12-30 00:00:00   | 360 d       |
 * | `-2`           | Gregorian, truncated  | 1968-05-23 00:00:00   | 365.25 d    |
 *                    Julian day (NASA)
 *
 * `time_ref` is passed in rather than read from a global, so these stay pure
 * functions; ROMS keeps the equivalent in `Rclock` in `mod_scalars`. The
 * reference date number is recomputed per call, which is a handful of integer
 * operations, rather than cached.
 *
 * The arithmetic is deliberately a transcription rather than a tidier
 * equivalent. Fortran `INT` and `AINT` truncate toward zero, as do a C++ cast
 * to `int` and `std::trunc`, and C++ integer division truncates toward zero
 * too, so each expression carries over term for term. Where ROMS uses `FLOOR`
 * -- and in the `time_ref = -2` conversion it does, which rounds differently
 * from `INT` on the negative side of the epoch -- `std::floor` is used here.
 *
 * Hours, minutes and seconds are not returned. ROMS derives them from the day
 * fraction through the tolerant rounding helper in `round_mod`, and
 * half-porting that would be worse than leaving it out; the day fraction is
 * available from remora_caldate's `yday` for a caller that needs it.
 *
 * Two quirks of the `time_ref = -2` calendar carry over, both ROMS's and both
 * confirmed against the Fortran rather than inferred:
 *
 *   - `datenum` splits the year as `y01 = MyYear/100`, which truncates toward
 *     zero, so it disagrees with the true Julian day number for years before 0.
 *     `dateclock.F` documents `datenum(-4713,11,24) = 0`; its own code returns
 *     1, and so does this. From year 1 on it agrees with the standard
 *     Fliegel-Van Flandern day number exactly.
 *   - `datevec` treats any day number below the 15 October 1582 changeover as a
 *     *truncated* one and adds the reference offset back, so it does not invert
 *     `datenum` before that date. That is the heuristic ROMS uses to accept both
 *     full and truncated Julian days in the same argument.
 *
 * Neither affects a model run in modern times, which is what this calendar is
 * for, but both would bite anyone reaching further back with it.
 */

//! Reference day number of the 15 October 1582 Gregorian changeover, as ROMS
//! spells it in datevec.
static constexpr int remora_gregorian_daynum = 2299161;

/** \brief Day of the year for a calendar date. ROMS `yearday`. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int remora_yearday (int year, int month, int day) noexcept
{
    const bool leap = (((year % 4) == 0) && ((year % 100) != 0)) || ((year % 400) == 0);
    const int fac = leap ? 1 : 2;
    return static_cast<int>((275.0 * month) / 9.0) - fac * ((month + 9) / 12) + day - 30;
}

/**
 * \brief Fractional day number for a calendar date. ROMS `datenum`.
 *
 * Reference values, per calendar:
 *   - proleptic Gregorian: `datenum(0000,01,01) = 1`, `datenum(0001,01,01) = 367`.
 *     The origin is Matlab's `datenum(0000,00,00) = 0`, hence the 61-day offset.
 *   - 360_day: `datenum(0000,01,01) = 0`, `datenum(0001,01,01) = 360`.
 *   - truncated Julian: `datenum(-4713,11,24) = 0`, `datenum(1968,05,23) = 2440000`.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real remora_datenum (amrex::Real time_ref, int year, int month, int day,
                            int hour = 0, int minutes = 0,
                            amrex::Real seconds = amrex::Real(0.0)) noexcept
{
    const int cal = static_cast<int>(time_ref);
    int my_day = 0;

    if (cal == -2) {
        // Julian day plus the Gregorian correction. Julian days formally start
        // at noon; here they start at midnight, so this runs 12 hours ahead of
        // the formal definition.
        int my_year;
        int my_month;
        if (month > 2) {
            my_year = year;
            my_month = month - 3;
        } else {
            my_year = year - 1;
            my_month = month + 9;
        }
        const int y01 = my_year / 100;
        my_year -= y01 * 100;
        my_day = (146097 * y01 / 4) + (1461 * my_year / 4) +
                 ((153 * my_month + 2) / 5) + day + 1721119;

    } else if (cal == -1) {
        // Every year is 360 days and every month 30.
        my_day = year * 360 + (month - 1) * 30 + (day - 1);

    } else {
        // Gregorian and proleptic Gregorian, year length 365.2425.
        constexpr int offset = 61;

        const int my_month = (month + 9) % 12;                          // Mar=0, ..., Feb=11
        const int my_year = year - static_cast<int>(0.1 * my_month);    // Jan or Feb: back one

        my_day = static_cast<int>(365.0 * my_year) +
                 static_cast<int>(0.25 * my_year) -
                 static_cast<int>(0.01 * my_year) +
                 static_cast<int>(0.0025 * my_year) +
                 static_cast<int>(0.1 * (my_month * 306.0 + 5.0)) +
                 (day - 1);

        if ((year == 0) && (month == 0) && (day == 0)) {
            my_day = 0;
        } else if (my_day < 0) {
            my_day += offset - 1;
        } else {
            my_day += offset;
        }
    }

    return amrex::Real(my_day) +
           amrex::Real(hour) / amrex::Real(24.0) +
           amrex::Real(minutes) / amrex::Real(1440.0) +
           seconds / amrex::Real(86400.0);
}

/**
 * \brief Day number of the reference date implied by \p time_ref.
 *
 * ROMS `ref_clock`, reduced to the one field the rest of this header needs:
 * `Rclock%DateNumber(1)`. The two negative calendars carry hardcoded values in
 * ROMS -- 359 for 360_day, whose epoch is shifted to 0000-12-30 to undo a
 * historical one-day offset, and 2440000 for the truncated Julian day.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real remora_ref_datenum (amrex::Real time_ref) noexcept
{
    const int cal = static_cast<int>(time_ref);

    if (cal > 0) {
        // Decode yyyymmdd.dd. The clamps are ROMS's, and they are what keeps a
        // malformed value from indexing off the end of the calendar.
        const int iyear = std::max(1, static_cast<int>(time_ref * amrex::Real(0.0001)));
        const int month = std::min(12, std::max(1, static_cast<int>(
            (time_ref - amrex::Real(iyear * 10000)) * amrex::Real(0.01))));
        const amrex::Real day = time_ref -
            std::trunc(time_ref * amrex::Real(0.01)) * amrex::Real(100.0);
        const int iday = std::max(1, static_cast<int>(day));
        const amrex::Real sec = (day - std::trunc(day)) * amrex::Real(86400.0);
        const int ihour = static_cast<int>(sec / amrex::Real(3600.0));
        const int minute = static_cast<int>(std::fmod(sec, amrex::Real(3600.0)) /
                                            amrex::Real(60.0));
        const int isec = static_cast<int>(std::fmod(sec, amrex::Real(60.0)));
        return remora_datenum(time_ref, iyear, month, iday, ihour, minute,
                              amrex::Real(isec));
    }
    if (cal == -1) {
        return amrex::Real(359.0);            // 0000-12-30 in the 360_day calendar
    }
    if (cal == -2) {
        return amrex::Real(2440000.0);        // 1968-05-23, truncated Julian offset
    }
    return remora_datenum(time_ref, 1, 1, 1); // 0001-01-01, equals 367
}

/**
 * \brief Calendar date for a day (or second) number. ROMS `datevec`.
 *
 * \p is_day_units selects fractional days over fractional seconds, as ROMS's
 * `IsDayUnits` does.
 *
 * ROMS's `time_ref = -2` branch can also run a proleptic *Julian* conversion,
 * but its `ProlepticJulian` switch is a local hardcoded to `.FALSE.`, so that
 * path is unreachable as shipped and is not ported.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void remora_datevec (amrex::Real time_ref, amrex::Real date_number, bool is_day_units,
                     int& year, int& month, int& day) noexcept
{
    const int cal = static_cast<int>(time_ref);

    if (cal == -2) {
        // A value at or past the Gregorian changeover is taken to be a full
        // Julian day number; anything smaller is a truncated one and gets the
        // reference offset added back.
        amrex::Real my_date_number;
        if (is_day_units) {
            my_date_number = (date_number >= amrex::Real(remora_gregorian_daynum))
                           ? date_number
                           : date_number + remora_ref_datenum(time_ref);
        } else {
            my_date_number =
                (date_number >= amrex::Real(remora_gregorian_daynum) * amrex::Real(86400.0))
                ? date_number / amrex::Real(86400.0)
                : date_number / amrex::Real(86400.0) + remora_ref_datenum(time_ref);
        }

        // Proleptic Gregorian, origin 24 November 4713 BC. FLOOR throughout,
        // not truncation: this calendar runs to dates before its own origin.
        amrex::Real jr = std::floor(my_date_number) - amrex::Real(1721119.0);
        amrex::Real js = amrex::Real(4.0) * jr - amrex::Real(1.0);
        amrex::Real yy = std::floor(js / amrex::Real(146097.0));
        jr = js - amrex::Real(146097.0) * yy;
        js = std::floor(jr * amrex::Real(0.25));
        js = amrex::Real(4.0) * js + amrex::Real(3.0);
        jr = std::floor(js / amrex::Real(1461.0));
        const amrex::Real dd =
            std::floor(((js - amrex::Real(1461.0) * jr) + amrex::Real(4.0)) * amrex::Real(0.25));
        js = amrex::Real(5.0) * dd - amrex::Real(3.0);
        const amrex::Real mo = std::floor(js / amrex::Real(153.0));
        yy = yy * amrex::Real(100.0) + jr;

        if (mo < amrex::Real(10.0)) {
            year = static_cast<int>(yy);
            month = static_cast<int>(mo + amrex::Real(3.0));
        } else {
            year = static_cast<int>(yy + amrex::Real(1.0));
            month = static_cast<int>(mo - amrex::Real(9.0));
        }
        day = static_cast<int>(((js - amrex::Real(153.0) * mo) + amrex::Real(5.0)) *
                               amrex::Real(0.2));

    } else if (cal == -1) {
        // 360_day: 12 months of 30 days, no leap years.
        int yday;
        if (is_day_units) {
            year = static_cast<int>(date_number / amrex::Real(360.0));
            yday = static_cast<int>(date_number - amrex::Real(year * 360) + amrex::Real(1.0));
        } else {
            year = static_cast<int>(date_number / amrex::Real(31104000.0));   // 360*86400
            yday = static_cast<int>((date_number - amrex::Real(year * 31104000) +
                                     amrex::Real(1.0)) / amrex::Real(86400.0));
        }
        month = ((yday - 1) / 30) + 1;
        day = ((yday - 1) % 30) + 1;

    } else {
        // Gregorian and proleptic Gregorian.
        constexpr amrex::Real offset = amrex::Real(61.0);

        amrex::Real my_date_number = is_day_units ? date_number
                                                  : date_number / amrex::Real(86400.0);
        if (my_date_number < offset) {          // adjust for the Matlab zero
            my_date_number = my_date_number - offset + amrex::Real(1.0);
        } else {                                // origin, datenum(0,0,0) = 0
            my_date_number = my_date_number - offset;
        }

        int my_year = static_cast<int>((10000.0 * std::trunc(my_date_number) + 14780.0) /
                                       3652425.0);
        auto days_before_year = [] (int y) {
            return static_cast<int>(365.0 * y) + static_cast<int>(0.25 * y) -
                   static_cast<int>(0.01 * y) + static_cast<int>(0.0025 * y);
        };
        int my_day = static_cast<int>(my_date_number) - days_before_year(my_year);
        if (my_day < 0) {                       // before 1 March; this keeps
            my_year -= 1;                       // leap years easy
            my_day = static_cast<int>(my_date_number) - days_before_year(my_year);
        }

        const int my_month = static_cast<int>((100.0 * my_day + 52.0) / 3060.0);
        month = ((my_month + 2) % 12) + 1;
        year = my_year + static_cast<int>((my_month + 2.0) / 12.0);
        day = my_day - static_cast<int>(0.1 * (my_month * 306.0 + 5.0)) + 1;

        // Matches Matlab "datestr" with the origin at 0000-00-00 00:00:00.
        if (date_number == amrex::Real(0.0)) {
            year = 0;
            month = 1;
            day = 0;
        }
    }
}

/**
 * \brief Model time in days to calendar date. ROMS `caldate`.
 *
 * \param[in]  time_ref           remora.time_ref; see the table at the top
 * \param[in]  current_time_days  model time in days, as ROMS `tdays`
 * \param[out] year               year including century
 * \param[out] month              month of the year, 1 = January
 * \param[out] day                day of the month
 * \param[out] yday               day of the year, including the day fraction,
 *                                so 00:00 on the first day of the year is
 *                                exactly 1.0
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void remora_caldate (amrex::Real time_ref, amrex::Real current_time_days,
                     int& year, int& month, int& day, amrex::Real& yday) noexcept
{
    const int cal = static_cast<int>(time_ref);
    const amrex::Real ref_date_number = remora_ref_datenum(time_ref);

    // For the truncated Julian day, a model time already past the reference is
    // taken to be an absolute day number rather than an elapsed interval.
    amrex::Real date_number;
    if ((cal == -2) && (current_time_days >= ref_date_number)) {
        date_number = current_time_days;
    } else {
        date_number = ref_date_number + current_time_days;
    }

    const amrex::Real day_fraction = std::abs(date_number - std::trunc(date_number));

    remora_datevec(time_ref, date_number, true, year, month, day);

    // The 360_day calendar has no leap years for yearday to reason about, so
    // ROMS takes its day of year straight from the day number.
    const int int_yday = (cal == -1)
        ? static_cast<int>(date_number - amrex::Real(year * 360) + amrex::Real(1.0))
        : remora_yearday(year, month, day);

    yday = amrex::Real(int_yday) + day_fraction;
}

/** \brief remora_caldate for callers that only want the year and day of year. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void remora_caldate (amrex::Real time_ref, amrex::Real current_time_days,
                     int& year, amrex::Real& yday) noexcept
{
    int month = 0;
    int day = 0;
    remora_caldate(time_ref, current_time_days, year, month, day, yday);
}

/**
 * \brief Whether \p time_ref names a calendar this header implements.
 *
 * ROMS's `CALENDAR` if-chain has no final `ELSE`, so a value below -2 leaves
 * the date undefined rather than failing. Callers should reject it up front.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool remora_time_ref_is_valid (amrex::Real time_ref) noexcept
{
    return static_cast<int>(time_ref) >= -2;
}

#endif

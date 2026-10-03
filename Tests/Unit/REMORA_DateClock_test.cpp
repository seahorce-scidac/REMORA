/**
 * \file REMORA_DateClock_test.cpp
 *
 * Unit tests for Source/Utils/REMORA_DateClock.H, REMORA's port of ROMS
 * `dateclock.F`.
 *
 * Every expectation comes from outside the implementation, so a failure means
 * the port is wrong rather than that it changed:
 *
 *   - Proleptic Gregorian day numbers are Matlab `datenum` values; ROMS
 *     documents its calendar as matching Matlab's.
 *   - The anchors -- `datenum(0001,01,01) = 367`, `datenum(1968,05,23) =
 *     2440000`, the 360_day pair -- are the values `dateclock.F` states in its
 *     own comments.
 *   - Reference dates and calendar names come off the four branches of ROMS
 *     `ref_clock`.
 *   - Days of the year are Gregorian arithmetic.
 *   - The round-trip checks are properties rather than tables: `datevec` must
 *     invert `datenum`, which needs no authority at all.
 *
 * The calendar works in double regardless of amrex::Real, so every case runs in
 * either build precision.
 */

#include "REMORA_DateClock.H"

#include <cstdio>
#include <string>

namespace {

int failures = 0;
int checks = 0;

void check_int (const char* what, int got, int expect)
{
    ++checks;
    if (got != expect) {
        ++failures;
        std::printf("FAIL %s: got %d, expected %d\n", what, got, expect);
    }
}

//! Exact comparison is deliberate: every value asserted below is an integer or
//! a short binary fraction, exactly representable in double, so the port should
//! reproduce it bit for bit rather than nearly.
void check_real (const char* what, double got, double expect)
{
    ++checks;
    if (got != expect) {
        ++failures;
        std::printf("FAIL %s: got %.6f, expected %.6f\n", what, got, expect);
    }
}

void check_str (const char* what, const std::string& got, const std::string& expect)
{
    ++checks;
    if (got != expect) {
        ++failures;
        std::printf("FAIL %s: got \"%s\", expected \"%s\"\n", what,
                    got.c_str(), expect.c_str());
    }
}

void check_date (const char* what, double time_ref, int year, int month, int day)
{
    int y = 0;
    int m = 0;
    int d = 0;
    remora_datevec(time_ref, remora_datenum(time_ref, year, month, day), true, y, m, d);
    ++checks;
    if ((y != year) || (m != month) || (d != day)) {
        ++failures;
        std::printf("FAIL %s: %04d-%02d-%02d round-tripped to %04d-%02d-%02d\n",
                    what, year, month, day, y, m, d);
    }
}

/** Day of the year. ROMS `yearday`. */
void test_yearday ()
{
    check_int("yearday(2020,01,01)", remora_yearday(2020, 1, 1), 1);
    check_int("yearday(2024,02,29)", remora_yearday(2024, 2, 29), 60);
    check_int("yearday(2024,03,01)", remora_yearday(2024, 3, 1), 61);
    check_int("yearday(2023,03,01)", remora_yearday(2023, 3, 1), 60);
    check_int("yearday(2024,12,31)", remora_yearday(2024, 12, 31), 366);
    check_int("yearday(2023,12,31)", remora_yearday(2023, 12, 31), 365);

    // 1900 is not a leap year and 2000 is: the century rules, not just the
    // four-year one.
    check_int("yearday(1900,03,01)", remora_yearday(1900, 3, 1), 60);
    check_int("yearday(2000,03,01)", remora_yearday(2000, 3, 1), 61);
}

/** Calendar date to day number. ROMS `datenum`. */
void test_datenum ()
{
    const double greg = 0.0;

    // dateclock.F's own reference values for this calendar.
    check_real("datenum(0000,01,01)", remora_datenum(greg, 0, 1, 1), 1.0);
    check_real("datenum(0001,01,01)", remora_datenum(greg, 1, 1, 1), 367.0);

    // Matlab datenum.
    check_real("datenum(2000,01,01)", remora_datenum(greg, 2000, 1, 1),
               730486.0);
    check_real("datenum(2020,01,01)", remora_datenum(greg, 2020, 1, 1),
               737791.0);
    check_real("datenum(2024,02,29)", remora_datenum(greg, 2024, 2, 29),
               739311.0);

    // The day fraction is hours/24 + minutes/1440 + seconds/86400.
    check_real("datenum(2020,01,01,12:00:00)",
               remora_datenum(greg, 2020, 1, 1, 12, 0, 0.0),
               737791.5);
    check_real("datenum(2020,01,01,06:00:00)",
               remora_datenum(greg, 2020, 1, 1, 6, 0, 0.0),
               737791.25);

    // 360_day: twelve 30-day months, no leap years. dateclock.F's values.
    const double day360 = -1.0;
    check_real("datenum(0000,01,01) [360_day]", remora_datenum(day360, 0, 1, 1),
               0.0);
    check_real("datenum(0001,01,01) [360_day]", remora_datenum(day360, 1, 1, 1),
               360.0);
    check_real("datenum(2020,01,01) [360_day]", remora_datenum(day360, 2020, 1, 1),
               727200.0);   // 2020 * 360

    // Truncated Julian day. dateclock.F documents datenum(1968,05,23) = 2440000.
    check_real("datenum(1968,05,23) [trunc Julian]",
               remora_datenum(-2.0, 1968, 5, 23), 2440000.0);
}

/** Day number of the reference date. ROMS `ref_clock`'s `Rclock%DateNumber(1)`. */
void test_ref_datenum ()
{
    // The three sentinel calendars, as dateclock.F hardcodes them.
    check_real("ref_datenum(0)", remora_ref_datenum(0.0),
               367.0);       // 0001-01-01
    check_real("ref_datenum(-1)", remora_ref_datenum(-1.0),
               359.0);       // 0000-12-30, 360_day
    check_real("ref_datenum(-2)", remora_ref_datenum(-2.0),
               2440000.0);   // 1968-05-23, truncated Julian

    // A yyyymmdd.dd date must decode to the same day number as naming its
    // fields outright.
    check_real("ref_datenum(20200101)", remora_ref_datenum(20200101.0),
               737791.0);
    check_real("ref_datenum(20200101.5)", remora_ref_datenum(20200101.5),
               737791.5);
    check_real("ref_datenum(20000101)", remora_ref_datenum(20000101.0),
               730486.0);
    check_real("ref_datenum(19680523)", remora_ref_datenum(19680523.0),
               718941.0);    // Matlab datenum(1968,5,23)

    // ROMS clamps a month or day of zero up to 1, so these are all 0001-01-01.
    check_real("ref_datenum(101)", remora_ref_datenum(101.0),
               367.0);
    check_real("ref_datenum(10101)", remora_ref_datenum(10101.0),
               367.0);
}

/** `datevec` must invert `datenum` over the range a model run can reach. */
void test_datevec_roundtrip ()
{
    // Proleptic Gregorian, not extended below the 1582 changeover: the header
    // documents that datevec reads a small day number as a truncated Julian one
    // when time_ref = -2, and that datenum(0,0,0) = 0 is a special case.
    const double greg = 0.0;
    for (int year = 1900; year <= 2100; year += 7) {
        for (int month = 1; month <= 12; ++month) {
            check_date("roundtrip [gregorian]", greg, year, month, 1);
            check_date("roundtrip [gregorian]", greg, year, month, 15);
            check_date("roundtrip [gregorian]", greg, year, month, 28);
        }
    }
    check_date("roundtrip 2024-02-29", greg, 2024, 2, 29);
    check_date("roundtrip 2000-02-29", greg, 2000, 2, 29);
    check_date("roundtrip 2024-12-31", greg, 2024, 12, 31);

    // 360_day, where every month has 30 days.
    const double day360 = -1.0;
    for (int year = 1900; year <= 2100; year += 7) {
        for (int month = 1; month <= 12; ++month) {
            check_date("roundtrip [360_day]", day360, year, month, 1);
            check_date("roundtrip [360_day]", day360, year, month, 30);
        }
    }

    // Truncated Julian, modern dates only.
    const double julian = -2.0;
    for (int year = 1900; year <= 2100; year += 7) {
        for (int month = 1; month <= 12; ++month) {
            check_date("roundtrip [trunc Julian]", julian, year, month, 1);
            check_date("roundtrip [trunc Julian]", julian, year, month, 28);
        }
    }
}

/** Model time to calendar date. ROMS `caldate`. */
void test_caldate ()
{
    int year = 0;
    int month = 0;
    int day = 0;
    double yday = 0.0;

    // Model time zero is the reference date itself: 0001-01-01 for time_ref = 0.
    remora_caldate(0.0, 0.0, year, month, day, yday);
    check_int("caldate(0, 0 d) year", year, 1);
    check_int("caldate(0, 0 d) month", month, 1);
    check_int("caldate(0, 0 d) day", day, 1);
    check_real("caldate(0, 0 d) yday", yday, 1.0);

    // 737424 days past 0001-01-01 is 1 January 2020 -- the number
    // remora.start_time = 63713433600 puts in the dstart output variable, so
    // this case ties that documented input to the calendar.
    remora_caldate(0.0, 737424.0, year, month, day, yday);
    check_int("caldate(0, 737424 d) year", year, 2020);
    check_int("caldate(0, 737424 d) month", month, 1);
    check_int("caldate(0, 737424 d) day", day, 1);
    check_real("caldate(0, 737424 d) yday", yday, 1.0);

    // Half a day in, the date is unchanged and yday carries the fraction.
    remora_caldate(0.0, 737424.5, year, month, day, yday);
    check_int("caldate(0, 737424.5 d) day", day, 1);
    check_real("caldate(0, 737424.5 d) yday", yday, 1.5);

    // A leap day reached from the reference date, not named directly:
    // datenum(2024,2,29) - datenum(0001,01,01) = 739311 - 367.
    remora_caldate(0.0, 738944.0, year, month, day, yday);
    check_int("caldate(0, 738944 d) year", year, 2024);
    check_int("caldate(0, 738944 d) month", month, 2);
    check_int("caldate(0, 738944 d) day", day, 29);
    check_real("caldate(0, 738944 d) yday", yday, 60.0);

    // 360_day. The epoch is 0000-12-30: ROMS shifts it back a day so day number
    // 0 lands on 01-Jan-0000, undoing a historical offset. So model time 1 day
    // is 0001-01-01 -- not 0, and not 360.
    remora_caldate(-1.0, 1.0, year, month, day, yday);
    check_int("caldate(-1, 1 d) year", year, 1);
    check_int("caldate(-1, 1 d) month", month, 1);
    check_int("caldate(-1, 1 d) day", day, 1);
    check_real("caldate(-1, 1 d) yday", yday, 1.0);

    // A year on. Every month has 30 days, so the last day of the year is the
    // 30th of the twelfth and yday is exactly 360.
    remora_caldate(-1.0, 360.0, year, month, day, yday);
    check_int("caldate(-1, 360 d) year", year, 1);
    check_int("caldate(-1, 360 d) month", month, 12);
    check_int("caldate(-1, 360 d) day", day, 30);
    check_real("caldate(-1, 360 d) yday", yday, 360.0);

    // The two-argument form must agree with the six-argument one.
    int year2 = 0;
    double yday2 = 0.0;
    remora_caldate(0.0, 737424.0, year2, yday2);
    check_int("caldate 2-arg year", year2, 2020);
    check_real("caldate 2-arg yday", yday2, 1.0);
}

/** The reference date behind `Rclock%string`. ROMS `ref_clock`. */
void test_ref_clock ()
{
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    int second = 0;

    // time_ref = 0: proleptic Gregorian from 0001-01-01 00:00:00.
    remora_ref_clock(0.0, year, month, day, hour, minute, second);
    check_int("ref_clock(0) year", year, 1);
    check_int("ref_clock(0) month", month, 1);
    check_int("ref_clock(0) day", day, 1);
    check_int("ref_clock(0) hour", hour, 0);
    check_int("ref_clock(0) minute", minute, 0);
    check_int("ref_clock(0) second", second, 0);

    // time_ref = -1: 360_day, epoch shifted back a day to 0000-12-30.
    remora_ref_clock(-1.0, year, month, day, hour, minute, second);
    check_int("ref_clock(-1) year", year, 0);
    check_int("ref_clock(-1) month", month, 12);
    check_int("ref_clock(-1) day", day, 30);

    // time_ref = -2: truncated Julian day from 1968-05-23.
    remora_ref_clock(-2.0, year, month, day, hour, minute, second);
    check_int("ref_clock(-2) year", year, 1968);
    check_int("ref_clock(-2) month", month, 5);
    check_int("ref_clock(-2) day", day, 23);

    remora_ref_clock(20200101.0, year, month, day, hour, minute, second);
    check_int("ref_clock(20200101) year", year, 2020);
    check_int("ref_clock(20200101) month", month, 1);
    check_int("ref_clock(20200101) day", day, 1);
    check_int("ref_clock(20200101) hour", hour, 0);

    // The fractional part is a time of day: .5 is 12:00, as dateclock.F's
    // "20020115.5 for 15 Jan 2002, 12:0:0" states.
    remora_ref_clock(20200101.5, year, month, day, hour, minute, second);
    check_int("ref_clock(20200101.5) day", day, 1);
    check_int("ref_clock(20200101.5) hour", hour, 12);
    check_int("ref_clock(20200101.5) minute", minute, 0);
    check_int("ref_clock(20200101.5) second", second, 0);

    remora_ref_clock(20250908.25, year, month, day, hour, minute, second);
    check_int("ref_clock(20250908.25) year", year, 2025);
    check_int("ref_clock(20250908.25) month", month, 9);
    check_int("ref_clock(20250908.25) day", day, 8);
    check_int("ref_clock(20250908.25) hour", hour, 6);
}

/** The CF time-units strings. ROMS `Rclock%string` and `Rclock%calendar`. */
void test_ref_strings ()
{
    check_str("ref_date_string(0)", remora_ref_date_string(0.0),
              "0001-01-01 00:00:00");
    check_str("ref_date_string(-1)", remora_ref_date_string(-1.0),
              "0000-12-30 00:00:00");
    check_str("ref_date_string(-2)", remora_ref_date_string(-2.0),
              "1968-05-23 00:00:00");

    // ref_clock assigns proleptic_gregorian to both non-negative branches.
    check_str("ref_calendar(0)", remora_ref_calendar(0.0),
              "proleptic_gregorian");
    check_str("ref_calendar(-1)", remora_ref_calendar(-1.0), "360_day");
    check_str("ref_calendar(-2)", remora_ref_calendar(-2.0), "gregorian");

    check_str("ref_date_string(20200101)",
              remora_ref_date_string(20200101.0), "2020-01-01 00:00:00");
    check_str("ref_date_string(20200101.5)",
              remora_ref_date_string(20200101.5), "2020-01-01 12:00:00");
    check_str("ref_date_string(20250908.25)",
              remora_ref_date_string(20250908.25), "2025-09-08 06:00:00");
    check_str("ref_calendar(20200101)",
              remora_ref_calendar(20200101.0), "proleptic_gregorian");
}

/** The guard on remora.time_ref. */
void test_time_ref_validation ()
{
    check_int("is_valid(0)", remora_time_ref_is_valid(0.0) ? 1 : 0, 1);
    check_int("is_valid(-1)", remora_time_ref_is_valid(-1.0) ? 1 : 0, 1);
    check_int("is_valid(-2)", remora_time_ref_is_valid(-2.0) ? 1 : 0, 1);
    check_int("is_valid(20200101)", remora_time_ref_is_valid(20200101.0) ? 1 : 0, 1);

    // ROMS's calendar if-chain has no final ELSE, so anything below -2 names no
    // calendar and has to be rejected up front.
    check_int("is_valid(-3)", remora_time_ref_is_valid(-3.0) ? 1 : 0, 0);
    check_int("is_valid(-100)", remora_time_ref_is_valid(-100.0) ? 1 : 0, 0);
}

} // namespace

int main ()
{
    test_yearday();
    test_datenum();
    test_ref_datenum();
    test_datevec_roundtrip();
    test_caldate();
    test_ref_clock();
    test_ref_strings();
    test_time_ref_validation();

    if (failures > 0) {
        std::printf("REMORA_DateClock: %d of %d checks FAILED\n", failures, checks);
        return 1;
    }

    std::printf("REMORA_DateClock: all %d checks passed\n", checks);
    return 0;
}

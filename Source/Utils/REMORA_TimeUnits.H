#ifndef REMORA_TimeUnits_H
#define REMORA_TimeUnits_H

#include <cctype>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "REMORA_DateClock.H"

/** \file REMORA_TimeUnits.H
 *
 * The CF `units` attribute of a NetCDF time variable, and how it maps onto the
 * model clock.
 *
 * A forcing, boundary, river or climatology file stamps its records with a
 * time variable whose `units` read `days since <date>`, or just `days`. The
 * model clock counts from the date `remora.time_ref` names (see
 * REMORA_DateClock.H), so a file time is on the model clock only after the
 * difference between the two reference dates is added:
 *
 *     model days = file days + (datenum(file reference) - datenum(time_ref))
 *
 * remora_time_units_offset_days computes that bracket, in the calendar
 * `time_ref` selects. It is zero when the file counts from the same date as the
 * run, which is the usual case, and a file with no reference date is taken to be
 * on the model clock already.
 *
 * remora_parse_time_units reads the attribute. It accepts the CF form,
 * `<unit> since YYYY-MM-DD[ hh:mm[:ss[.s]]][ UTC|Z|+00:00]`, and it skips words
 * between `since` and the date, so `days since initialization, 2020-01-01
 * 00:00:00` reads as 2020-01-01. It does not interpret a time zone other than
 * UTC. Host only: it works on std::string.
 */

//! A parsed CF time-units attribute.
struct REMORATimeUnits
{
    //! The unit word, lower-cased: `days`, `seconds`, ...
    std::string unit;
    //! Whether a `since <date>` reference date followed the unit
    bool has_reference = false;
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    double second = 0.0;
};

namespace remora_time_units_detail {

inline std::string lower (std::string s)
{
    for (auto& c : s) {
        // A NUL terminator stored with the attribute text would otherwise glue
        // itself to the last word.
        c = (c == '\0') ? ' ' : static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

//! `hh:mm[:ss[.s]]`, optionally suffixed `z`. Fills the time fields on success.
inline bool parse_clock (const std::string& tok, REMORATimeUnits& out)
{
    int h = 0, m = 0, n = 0;
    double s = 0.0;
    bool matched = (std::sscanf(tok.c_str(), "%d:%d:%lf%n", &h, &m, &s, &n) == 3);
    if (!matched) {
        s = 0.0;
        n = 0;
        matched = (std::sscanf(tok.c_str(), "%d:%d%n", &h, &m, &n) == 2);
    }
    if (!matched) {
        return false;
    }
    const std::string rest = tok.substr(n);
    if (!rest.empty() && rest != "z") {
        return false;
    }
    out.hour = h;
    out.minute = m;
    out.second = s;
    return true;
}

//! `[-]YYYY-MM-DD`, optionally followed by `T` and a clock. Fills the date
//! fields on success; a clock attached with `T` is parsed too.
inline bool parse_date (const std::string& tok, REMORATimeUnits& out)
{
    int y = 0, mo = 0, d = 0, n = 0;
    if (std::sscanf(tok.c_str(), "%d-%d-%d%n", &y, &mo, &d, &n) != 3) {
        return false;
    }
    const std::string rest = tok.substr(n);
    if (!rest.empty()) {
        REMORATimeUnits clock;
        if (rest[0] != 't' || !parse_clock(rest.substr(1), clock)) {
            return false;
        }
        out.hour = clock.hour;
        out.minute = clock.minute;
        out.second = clock.second;
    }
    out.year = y;
    out.month = mo;
    out.day = d;
    return true;
}

//! A time-zone designator that means UTC
inline bool is_utc (const std::string& tok)
{
    static const char* const names[] = {"utc", "gmt", "z", "+0", "-0", "+00", "-00",
                                        "+0:00", "-0:00", "+00:00", "-00:00",
                                        "+0000", "-0000"};
    for (const char* name : names) {
        if (tok == name) { return true; }
    }
    return false;
}

} // namespace remora_time_units_detail

/**
 * \brief Parse a CF time-units attribute.
 *
 * @param[in]  units  the attribute text, e.g. `days since 2020-01-01 00:00:00`
 * @param[out] out    the unit word and, if given, the reference date
 * @param[out] error  why parsing failed, when it did
 * @return whether the attribute parsed. The unit word is not checked here; the
 *         caller decides which units it accepts.
 */
inline bool remora_parse_time_units (const std::string& units, REMORATimeUnits& out,
                                     std::string& error)
{
    using namespace remora_time_units_detail;

    out = REMORATimeUnits{};
    error.clear();

    std::vector<std::string> tok;
    {
        std::istringstream iss(lower(units));
        std::string w;
        while (iss >> w) { tok.push_back(w); }
    }
    if (tok.empty()) {
        error = "the units attribute is empty";
        return false;
    }
    out.unit = tok[0];

    std::size_t i = 1;
    while (i < tok.size() && tok[i] != "since") { ++i; }
    if (i == tok.size()) {
        return true;                       // a bare unit: no reference date
    }

    // The first word after "since" that reads as a date. Anything before it, such
    // as "initialization,", is commentary.
    for (++i; i < tok.size(); ++i) {
        if (parse_date(tok[i], out)) { break; }
    }
    if (i == tok.size()) {
        error = "no date of the form yyyy-mm-dd follows 'since'";
        return false;
    }
    out.has_reference = true;
    ++i;

    if (i < tok.size() && parse_clock(tok[i], out)) { ++i; }
    if (i < tok.size()) {
        if (!is_utc(tok[i])) {
            error = "'" + tok[i] + "' is not a time of day or UTC; reference dates in"
                    " another time zone are not supported";
            return false;
        }
        ++i;
    }
    if (i < tok.size()) {
        error = "unexpected '" + tok[i] + "' after the reference date";
        return false;
    }

    if (out.month < 1 || out.month > 12 || out.day < 1 || out.day > 31 ||
        out.hour < 0 || out.hour > 23 || out.minute < 0 || out.minute > 59 ||
        out.second < 0.0 || out.second >= 61.0) {
        error = "the reference date is out of range";
        return false;
    }
    return true;
}

/**
 * \brief Days to add to a time counted from \p units' reference date to count it
 *        from the date \p time_ref names instead.
 *
 * Zero for a file with no reference date, which is taken to be on the model
 * clock already, and zero when the two dates agree. The day numbers are both
 * taken in the calendar \p time_ref selects, as remora_caldate does.
 */
inline double remora_time_units_offset_days (const REMORATimeUnits& units, double time_ref)
{
    if (!units.has_reference) {
        return 0.0;
    }
    return remora_datenum(time_ref, units.year, units.month, units.day,
                          units.hour, units.minute, units.second)
         - remora_ref_datenum(time_ref);
}

/**
 * \brief Whether a CF `calendar` attribute names the calendar \p time_ref selects.
 *
 * `360_day` goes with `time_ref = -1` and nothing else. Every other `time_ref`
 * is Gregorian: `gregorian`, `standard` and `proleptic_gregorian` all agree for
 * dates after 1582, which is where the model clock lives. `julian`, `noleap`
 * and the rest count days differently and are rejected.
 */
inline bool remora_calendar_matches_time_ref (const std::string& calendar, double time_ref)
{
    const std::string c = remora_time_units_detail::lower(calendar);
    // The attribute may carry a NUL terminator, which lower() turned into a space.
    const std::string name = c.substr(0, c.find_last_not_of(' ') + 1);
    if (static_cast<int>(time_ref) == -1) {
        return name == "360_day";
    }
    return name == "gregorian" || name == "standard" || name == "proleptic_gregorian";
}

/** \brief The reference date of \p units as `YYYY-MM-DD hh:mm:ss`. */
inline std::string remora_time_units_date_string (const REMORATimeUnits& units)
{
    char buf[6 * 11 + 5 + 1];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d",
                  units.year, units.month, units.day, units.hour, units.minute,
                  static_cast<int>(units.second));
    return std::string(buf);
}

#endif

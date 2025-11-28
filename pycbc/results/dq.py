'''This module contains utilities for following up search triggers'''

# JavaScript for searching the aLOG
redirect_javascript = """<script type="text/javascript">
function redirect(form,way)
{
        // Set location to form and submit.
        if(form != '')
        {
                document.forms[form].action=way;
                document.forms[form].submit();
        }
        else
        {
                window.top.location = way;
        }
}
</script>"""

search_form_string="""<form name="%s_alog_search" id="%s_alog_search" method="post">
<input type="hidden" name="srcDateFrom" id="srcDateFrom" value="%s" size="20"/>
<input type="hidden" name="srcDateTo" id="srcDateTo" value="%s" size="20"/>
</form>"""

data_h1_string = """
<a href=https://ldas-jobs.ligo-wa.caltech.edu/~detchar/summary/day/%s>
Summary</a>
&nbsp;
<a onclick="redirect('h1_alog_search',
'https://alog.ligo-wa.caltech.edu/aLOG/includes/search.php?adminType=search');
return true;">aLOG</a>"""

data_l1_string="""
<a href=https://ldas-jobs.ligo-la.caltech.edu/~detchar/summary/day/%s>
Summary</a>
&nbsp;
<a onclick="redirect('l1_alog_search',
'https://alog.ligo-la.caltech.edu/aLOG/includes/search.php?adminType=search');
return true;">aLOG</a>"""


def get_summary_page_link(ifo, utc_time):
    """Return a string that links to the summary page and aLOG for this ifo

    Parameters
    ----------
    ifo : string
        The detector name
    utc_time : datetime.date or datetime.datetime or sequence
        Either a datetime.date/datetime.datetime object, or a sequence whose
        first three elements are year, month, day (in that order).

    Returns
    -------
    return_string : string
        String containing HTML for links to summary page and aLOG search
    """
    search_form = search_form_string
    data = {'H1': data_h1_string, 'L1': data_l1_string}
    if ifo not in data:
        return ifo
    else:
        # support datetime/date objects or sequences (year, month, day)
        try:
            if hasattr(utc_time, 'year') and hasattr(utc_time, 'month') and hasattr(utc_time, 'day'):
                year = int(utc_time.year)
                month = int(utc_time.month)
                day = int(utc_time.day)
            else:
                year = int(utc_time[0])
                month = int(utc_time[1])
                day = int(utc_time[2])
        except Exception:
            raise TypeError("utc_time must be a datetime/date or a sequence (year, month, day)")

        # alog format is day-month-year
        alog_utc = '%02d-%02d-%4d' % (day, month, year)
        # summary page is exactly the reverse
        ext = '%4d%02d%02d' % (year, month, day)
        return_string = search_form % (ifo.lower(), ifo.lower(), alog_utc, alog_utc)
        return return_string + data[ifo] % ext


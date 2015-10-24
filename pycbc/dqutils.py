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

def get_summary_page_link(ifo, utc_time):
    """Return a string that links to the summary page and aLOG for this ifo

    Parameters
    ----------
    ifo : string
        The detector name
    utc_time : sequence
        First three elements must be strings giving year, month, day resp.

    Returns
    -------
    return_string : string
        String containing HTML for links to summary page and aLOG search
    """
    search_form = \
        """<form name="%s_alog_search" id="%s_alog_search" method="post">
<input type="hidden" name="srcDateFrom" id="srcDateFrom" value="%s" size="20"/>
<input type="hidden" name="srcDateTo" id="srcDateTo" value="%s" size="20"/>
           </form>"""
    data = {'H1':"""H1&nbsp;<a href=https://ldas-jobs.ligo-wa.caltech.edu/~detchar/summary/day/%s>Summary</a>&nbsp;<a onclick="redirect('h1_alog_search','https://alog.ligo-wa.caltech.edu/aLOG/includes/search.php?adminType=search'); return true;">aLOG</a>""",
            'L1':"""L1&nbsp;<a href=https://ldas-jobs.ligo-la.caltech.edu/~detchar/summary/day/%s>Summary</a>&nbsp;<a onclick="redirect('l1_alog_search','https://alog.ligo-la.caltech.edu/aLOG/includes/search.php?adminType=search'); return true;">aLOG</a>""" }
    if ifo not in data:
        return ifo
    else:
        # alog format is day-month-year
        alog_utc = '%02d-%02d-%4d' % (utc_time[2], utc_time[1], utc_time[0])
        # summary page is exactly the reverse
        ext = '%4d%02d%02d' % (utc_time[0], utc_time[1], utc_time[2])
        return_string = search_form % (ifo.lower(), ifo.lower(), alog_utc, alog_utc)
        return return_string + data[ifo] % ext


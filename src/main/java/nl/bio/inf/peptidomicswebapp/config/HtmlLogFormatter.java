package nl.bio.inf.peptidomicswebapp.config;

import java.util.Date;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.LogRecord;

/**
 *  This class setups the table for the HTML Log file.
 * @author Jan Alfonso Busker
 */

public class HtmlLogFormatter extends Formatter {
    /**
     *
     * @param record the log record to be formatted.
     * @return a table row with the information
     */
    public String format(LogRecord record) {
        return ("<tr><td>"
                + (new Date(record.getMillis()))
                + "</td><td>"
                + record.getLevel()
                + "</td><td>"
                + record.getSourceClassName()
                + "</td><td>"
                + record.getMessage()
                + "</td></tr>\n");
    }

    /**
     *
     * @param h  The target handler (can be null)
     * @return String containing the header for the table
     */
    public String getHead(Handler h) {
        return ("""
                <html>
                <head> <meta http-equiv='Content-Type' content='text/html; charset=UTF-8'/> </head><body>
                <Table border>
                <tr><td>Time</td><td>Log Level</td><td>Class Name</td><td>Log Message</td></tr>
                """);
    }

    /**
     *
     * @param h  The target handler (can be null)
     * @return String containing the tail/end for the table
     */
    public String getTail(Handler h) {
        return ("</table>\n</body>\n</html>");
    }
}

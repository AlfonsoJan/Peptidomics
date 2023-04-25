package nl.bio.inf.peptidomicswebapp.config;

import java.util.Date;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.LogRecord;

/**
 *  This class setups the table for the HTML Log file.
 * @Jan Alfonso
 */

public class HtmlLogFormatter extends Formatter {
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

    public String getHead(Handler h) {
        return ("""
                <html>
                <head> <meta http-equiv='Content-Type' content='text/html; charset=UTF-8'/> </head><body>
                <Table border>
                <tr><td>Time</td><td>Log Level</td><td>Class Name</td><td>Log Message</td></tr>
                """);
    }

    public String getTail(Handler h) {
        return ("</table>\n</body>\n</html>");
    }
}

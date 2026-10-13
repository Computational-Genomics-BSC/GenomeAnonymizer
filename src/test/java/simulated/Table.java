package simulated;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Minimal tab separated table with a header, used for the truth files, the validation reports and the baselines.
 * The first column is the key of each row
 * @author Nicolas Gaitan
 */
public class Table {

    public static final String FIELD_SEPARATOR = "\t";
    public static final String HEADER_PREFIX = "#";

    private final List<String> header;
    private final Map<String, Map<String, String>> rows = new LinkedHashMap<>();

    public Table(String... header) {
        this.header = new ArrayList<>(Arrays.asList(header));
    }

    public List<String> getHeader() {
        return header;
    }

    public void addRow(Object... values) {
        if (values.length != header.size()) {
            throw new IllegalArgumentException("Row with " + values.length + " values for a table with " + header.size() + " columns");
        }
        Map<String, String> row = new LinkedHashMap<>();
        for (int i = 0; i < values.length; i++) {
            row.put(header.get(i), String.valueOf(values[i]));
        }
        rows.put(String.valueOf(values[0]), row);
    }

    public Map<String, Map<String, String>> getRows() {
        return rows;
    }

    public Map<String, String> getRow(String key) {
        return rows.get(key);
    }

    public int size() {
        return rows.size();
    }

    public void write(File file) throws IOException {
        file.getAbsoluteFile().getParentFile().mkdirs();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(HEADER_PREFIX + String.join(FIELD_SEPARATOR, header));
            for (Map<String, String> row : rows.values()) {
                writer.println(String.join(FIELD_SEPARATOR, row.values()));
            }
        }
    }

    public static Table read(File file) throws IOException {
        Table answer = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line = reader.readLine();
            while (line != null) {
                if (answer == null) {
                    if (!line.startsWith(HEADER_PREFIX)) throw new IOException("Missing header in " + file);
                    answer = new Table(line.substring(HEADER_PREFIX.length()).split(FIELD_SEPARATOR, -1));
                }
                else if (!line.isEmpty()) {
                    answer.addRow((Object[]) line.split(FIELD_SEPARATOR, -1));
                }
                line = reader.readLine();
            }
        }
        if (answer == null) throw new IOException("Empty table file " + file);
        return answer;
    }
}

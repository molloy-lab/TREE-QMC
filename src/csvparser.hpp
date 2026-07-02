#ifndef CSVPARSER_HPP
#define CSVPARSER_HPP
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <variant>
#include <algorithm>
#include <iomanip>
using std::cout;
using Cell = std::variant<int, float, std::string>;
using Row = std::vector<Cell>;

inline constexpr char delimiters[] = {',', ';', '|', '\t'};
enum class FieldType
{
    INTEGER,
    STRING,
    FLOAT
};
struct fieldstruct
{
    FieldType key;
    std::string value;
    bool missingdata = false;
};

inline FieldType Ftypeidentifier(const std::string &value)
{
    if (value.empty())
    {
        cout << "\033[33mWarning\033[0m : Might have missing value in 1st data row of CSV. This will cause that field's datatype to be inferred as STRING.\n";
        return FieldType::STRING;
    }

    try
    {
        size_t idx;
        int i = std::stoi(value, &idx);
        if (idx == value.length())
            return FieldType::INTEGER; // Entire string parsed
    }
    catch (...)
    {
    }

    try
    {
        size_t idx;
        float f = std::stof(value, &idx);
        if (idx == value.length())
            return FieldType::FLOAT; // Entire string parsed
    }
    catch (...)
    {
    }

    return FieldType::STRING;
}
inline char delimiteridentifier(const std::string &line)
{
    for (char delimiter : delimiters)
    {
        if (line.find(delimiter) != std::string::npos)
        {
            return delimiter;
        }
    }
    cout << "\033[31mError\033[0m : No valid delimiter found in the file.\n";
    return ',';
}

class DataFrame
{
private:
    std::vector<fieldstruct> fieldsidentifier(std::ifstream &file, char delimiter=',')
    {
        std::vector<fieldstruct> fields;
        std::string line;
        bool inside_quotes = false;
        int pos = 0, prevpos = 0;
        fieldstruct field;
        std::string value;
        if (!file.is_open())
        {
            cout << "File not found\n";
        }
        else
        {
            std::getline(file, line);
            std::getline(file, line);
            line += delimiter;
            for (char ch : line)
            {
                if (ch == '"')
                    inside_quotes = !inside_quotes;
                if (ch == delimiter && !inside_quotes)
                {
                    value.clear();
                    fieldcount++;
                    for (int i = prevpos; i < pos; i++)
                    {
                        value += line[i];
                    }
                    field.key = Ftypeidentifier(value);
                    field.value = value;
                    fields.push_back(field);
                    prevpos = pos + 1;
                }
                pos++;
            }
            // move the file pointer to the start of the file
            file.clear();
            file.seekg(0, std::ios::beg);
            std::getline(file, line);
            line += delimiter;
            pos = 0;
            prevpos = 0;
            int j = 0;
            for (char ch : line)
            {
                if (ch == '"')
                    inside_quotes = !inside_quotes;
                if (ch == delimiter && !inside_quotes)
                {
                    value.clear();
                    for (int i = prevpos; i < pos; i++)
                    {
                        value += line[i];
                    }
                    fields[j].value = value;
                    prevpos = pos + 1;
                    j++;
                }
                pos++;
            }
        }
        if (fields.size() == 1)
        {
            cout << "\n\033[0;33mWarning\033[0m : Only one field is identified.\n";
            cout << "Check the delimiter used in the file and modify your command line according to it\n";
        }
        return fields;
    }

    void datacounter(std::ifstream &file)
    {
        std::string line;
        rowscount = 0;
        do
        {
            line.clear();
            std::getline(file, line);
            rowscount++;
        } while (!line.empty());
        rowscount--;
    }

    void createDataFrame(std::ifstream &file, char delimiter=',')
    {
        Row row;
        Cell cell;
        std::string line;
        std::string value;
        bool inside_quotes = false;
        std::getline(file, line); // Skip the header line
        line.clear();
        while (std::getline(file, line))
        {
            int pos = 0, prevpos = 0;
            int column = 0;
            line += delimiter;
            for (char ch : line)
            {
                if (ch == '"')
                    inside_quotes = !inside_quotes;
                if (ch == delimiter && !inside_quotes)
                {
                    value.clear();
                    for (int i = prevpos; i < pos; i++)
                    {
                        value += line[i];
                    }
                    if (TrueFields[column].key == FieldType::INTEGER)
                    {
                        if (value.empty())
                        {
                            TrueFields[column].missingdata = true;
                            cell = 0;
                        }
                        else
                        {
                            try
                            {
                                cell = std::stoi(value);
                            }
                            catch (const std::invalid_argument &)
                            {
                                cell = 0;
                            }
                        }
                    }
                    else if (TrueFields[column].key == FieldType::FLOAT)
                    {
                        if (value.empty())
                        {
                            TrueFields[column].missingdata = true;
                            cell = 0.0f;
                        }
                        else
                        {
                            try
                            {
                                cell = std::stof(value);
                            }
                            catch (const std::invalid_argument &)
                            {
                                cell = 0.0f;
                            }
                        }
                    }
                    else
                    {
                        if (value.empty())
                        {
                            TrueFields[column].missingdata = true;
                            cell = "NA";
                        }
                        else
                        {
                            cell = value;
                        }
                    }
                    row.push_back(cell);
                    column++;
                    prevpos = pos + 1;
                }
                pos++;
            }
            data.push_back(row);
            row.clear();
            line.clear();
        }
    }
  

public:
    std::vector<Row> data;
    std::vector<fieldstruct> TrueFields;
    int fieldcount = 0;
    int rowscount = 0;

    DataFrame() {}

    void reset()
    {
        data.clear();
        TrueFields.clear();
        fieldcount = 0;
        rowscount = 0;
    }

    void read_csv(const std::string &filename, char delimiter=',')
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            cout << "Error opening file.\n";
            return;
        }
        else
        {
            std::string line;
            std::getline(file, line);
            delimiter = delimiteridentifier(line);
            file.clear();
            file.seekg(0, std::ios::beg);
        }
        TrueFields = fieldsidentifier(file, delimiter);
        datacounter(file);
        file.clear();
        file.seekg(0, std::ios::beg);
        fieldcount = TrueFields.size();
        createDataFrame(file, delimiter);
    }

    void save(const std::string &filename)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            cout << "Error opening file for writing.\n";
            return;
        }
        for (int i = 0; i < fieldcount; ++i)
        {
            file << TrueFields[i].value;
            if (i < fieldcount - 1)
            {
                file << ",";
            }
        }
        file << "\n";
        for (const auto &row : data)
        {
            for (int i = 0; i < fieldcount; ++i)
            {
                if (std::holds_alternative<int>(row[i]))
                {
                    file << std::get<int>(row[i]);
                }
                else if (std::holds_alternative<float>(row[i]))
                {
                    file << std::get<float>(row[i]);
                }
                else
                {
                    file << std::get<std::string>(row[i]);
                }
                if (i < fieldcount - 1)
                {
                    file << ",";
                }
            }
            file << "\n";
        }
        file.close();
    }

    void copy(DataFrame &new_df,const std::vector<std::string> &field_names){
        std::vector<int> field_indices;
        if (field_names.empty())
        {
            cout << "\033[31mError\033[0m : No fields specified for the new DataFrame.\n";
            return;
        }
        bool flag;
        for(const std::string &name : field_names){
            flag = false;
            for(int i = 0; i < TrueFields.size(); ++i){
                if(TrueFields[i].value == name){
                    field_indices.push_back(i);
                    flag = true;
                    break;
                }
            }
            if(!flag){
                cout << "\033[31mError\033[0m : Field '" << name << "' not found in the DataFrame.\n";
                return;
            }
        }
        new_df.fieldcount = field_indices.size();
        new_df.TrueFields.resize(new_df.fieldcount);
        for(int i = 0; i < new_df.fieldcount; ++i){
            new_df.TrueFields[i] = TrueFields[field_indices[i]];
        }
        new_df.rowscount = rowscount;
        new_df.data.resize(rowscount);
        for(int i = 0; i < rowscount; ++i){
            new_df.data[i].resize(new_df.fieldcount);
            for(int j = 0; j < new_df.fieldcount; ++j){
                new_df.data[i][j] = data[i][field_indices[j]];
            }
        }
    }

    void copy(DataFrame &new_df, const std::vector<int> &field_indices)
    {
        if (field_indices.empty())
        {
            cout << "\033[31mError\033[0m : No fields specified for the new DataFrame.\n";
            return;
        }
        new_df.fieldcount = field_indices.size();
        new_df.TrueFields.resize(new_df.fieldcount);
        for (int i = 0; i < new_df.fieldcount; ++i)
        {
            if (field_indices[i] < 0 || field_indices[i] >= TrueFields.size())
            {
                cout << "\033[31mError\033[0m : Invalid field index " << field_indices[i] << ".\n";
                return;
            }
            new_df.TrueFields[i] = TrueFields[field_indices[i]];
        }
        new_df.rowscount = rowscount;
        new_df.data.resize(rowscount);
        for (int i = 0; i < rowscount; ++i)
        {
            new_df.data[i].resize(new_df.fieldcount);
            for (int j = 0; j < new_df.fieldcount; ++j)
            {
                new_df.data[i][j] = data[i][field_indices[j]];
            }
        }
    }

    void add_row(const Row &new_row)
    {
        if (new_row.size() != fieldcount)
        {
            cout << "\033[31mError\033[0m : Row size does not match the number of fields.\n";
            return;
        }
        data.push_back(new_row);
        rowscount++;
    }

    void change_value(std::string field_name, int row_index, const Cell &new_value)
    {
        if (row_index < 0 || row_index >= rowscount)
        {
            cout << "Invalid row index.\n";
            return;
        }
        auto it = std::find_if(TrueFields.begin(), TrueFields.end(),
                               [&field_name](const fieldstruct &field)
                               { return field.value == field_name; });
        if (it == TrueFields.end())
        {
            cout << "Field not found.\n";
            return;
        }
        int column_index = std::distance(TrueFields.begin(), it);
        if (column_index < 0 || column_index >= fieldcount)
        {
            cout << "Invalid column index.\n";
            return;
        }
        data[row_index][column_index] = new_value;
    }

    template <typename T>
    void get_column(std::vector<T> &x, const int index)
    {
        if (index < 0 || index >= fieldcount)
        {
            cout << "Invalid column index.\n";
        }
        for (int i = 0; i < rowscount; ++i)
        {
            if (std::holds_alternative<T>(data[i][index]))
            {
                x.push_back(std::get<T>(data[i][index]));
            }
            else
            {
                cout << "Type mismatch in column " << index << ".\n";
            }
        }
    }

    void sort(int column_index, bool ascending=true)
    {
        if (column_index < 0 || column_index >= fieldcount)
        {
            cout << "Invalid column index.\n";
            return;
        }
        std::sort(data.begin(), data.end(), [column_index, ascending](const Row &a, const Row &b)
                  {
            const Cell& cell_a = a[column_index];
            const Cell& cell_b = b[column_index];
            if (std::holds_alternative<int>(cell_a) && std::holds_alternative<int>(cell_b)) {
                return ascending ? std::get<int>(cell_a) < std::get<int>(cell_b)
                                : std::get<int>(cell_a) > std::get<int>(cell_b);
            } else if (std::holds_alternative<float>(cell_a) && std::holds_alternative<float>(cell_b)) {
                return ascending ? std::get<float>(cell_a) < std::get<float>(cell_b)
                                : std::get<float>(cell_a) > std::get<float>(cell_b);
            } else {
                return ascending ? std::get<std::string>(cell_a) < std::get<std::string>(cell_b)
                                : std::get<std::string>(cell_a) > std::get<std::string>(cell_b);
            } });
    }

    void print()
    {
        std::string value;
        FieldType key;
        int maxlen[fieldcount] = {0};
        for (int i = 0; i < fieldcount; i++)
        {
            key = TrueFields[i].key;
            if (key == FieldType::INTEGER)
            {
                for (int j = 0; j < rowscount; j++)
                {
                    value = std::to_string(data[j][i].index() == 0 ? std::get<int>(data[j][i]) : 0);
                    if (value.length() > maxlen[i])
                        maxlen[i] = value.length();
                }
            }
            else if (key == FieldType::FLOAT)
            {
                for (int j = 0; j < rowscount; j++)
                {
                    value = std::to_string(data[j][i].index() == 1 ? std::get<float>(data[j][i]) : 0.0f);
                    value.erase(value.find_last_not_of('0') + 1, std::string::npos);
                    if (!value.empty() && value.back() == '.')
                        value.pop_back();
                    if (value.length() > maxlen[i])
                        maxlen[i] = value.length();
                }
            }
            else
            {
                for (int j = 0; j < rowscount; j++)
                {
                    value = data[j][i].index() == 2 ? std::get<std::string>(data[j][i]) : "";
                    if (value.length() > maxlen[i])
                        maxlen[i] = value.length();
                }
            }
        }
        for (int i = 0; i < fieldcount; i++)
        {
            if (TrueFields[i].value.length() > maxlen[i])
            {
                maxlen[i] = TrueFields[i].value.length();
            }
        }
        cout << "\n+";
        for (int i = 0; i < fieldcount; i++)
        {
            cout << std::string(maxlen[i] + 2, '-') << "+";
        }
        cout << "\n";
        cout << "|";
        for (int i = 0; i < fieldcount; i++)
        {
            cout << " " << std::left << std::setw(maxlen[i]) << TrueFields[i].value << " |";
        }
        cout << "\n";
        cout << "+";
        for (int i = 0; i < fieldcount; i++)
        {
            cout << std::string(maxlen[i] + 2, '-') << "+";
        }
        cout << "\n";
        for (int j = 0; j < rowscount; j++)
        {
            cout << "|";
            for (int i = 0; i < fieldcount; i++)
            {
                key = TrueFields[i].key;
                if (key == FieldType::INTEGER)
                {
                    value = std::to_string(data[j][i].index() == 0 ? std::get<int>(data[j][i]) : 0);
                }
                else if (key == FieldType::FLOAT)
                {
                    value = std::to_string(data[j][i].index() == 1 ? std::get<float>(data[j][i]) : 0.0f);
                    value.erase(value.find_last_not_of('0') + 1, std::string::npos);
                    if (!value.empty() && value.back() == '.')
                        value.pop_back();
                }
                else
                {
                    value = data[j][i].index() == 2 ? std::get<std::string>(data[j][i]) : "";
                }
                cout << " " << std::left << std::setw(maxlen[i]) << value << " |";
            }
            cout << "\n";
        }
        cout << "+";
        for (int i = 0; i < fieldcount; i++)
        {
            cout << std::string(maxlen[i] + 2, '-') << "+";
        }
        cout << "\n";
    }

    void print_info()
    {
        cout << "DataFrame Information:\n";
        cout << "Number of Fields: " << fieldcount << "\n";
        cout << "Number of Rows: " << rowscount << "\n";
        cout << "\nField Information:\n";
        for (int i = 0; i < fieldcount; i++)
        {
            cout << "Field " << i + 1 << ": " << TrueFields[i].value
                 << " (Type: " << (TrueFields[i].key == FieldType::INTEGER ? "INTEGER" : TrueFields[i].key == FieldType::FLOAT ? "FLOAT"
                                                                                                         : "STRING")
                 << ", Missing Data: " << (TrueFields[i].missingdata ? "Yes" : "No") << ")\n";
        }
    }
};

#endif
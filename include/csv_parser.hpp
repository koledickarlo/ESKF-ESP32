#pragma once

#include <fstream>
#include <iostream>
#include <iterator>
#include <measurement_types.hpp>
#include <sstream>
#include <string>
#include <vector>

class CsvParser {
  public:
    inline std::string operator[](std::size_t index) const {
        return std::string(&m_line[m_data[index] + 1], m_data[index + 1] - (m_data[index] + 1));
    }
    inline std::size_t size() const { return m_data.size() - 1; }
    void readNextRow(std::istream &str);

  private:
    std::string m_line;
    std::vector<int> m_data;
};

std::istream &operator>>(std::istream &str, CsvParser &data);

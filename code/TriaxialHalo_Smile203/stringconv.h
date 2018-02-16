#pragma once
#include <string>
#include <sstream>
#include <iomanip>

/**
* The StringVariant class is a simple string-based variant implementation that allows
* the user to easily convert between simple numeric/string types.
*/
namespace smile{

class StringVariant
{
    std::string data;
public:
    StringVariant() : data() {}
    StringVariant(const std::string &src) : data(src) {};

    template<typename ValueType>
    StringVariant(ValueType val)
    {
        std::ostringstream stream;
        stream << val;
        data.assign(stream.str());
    }

    template<typename ValueType>
    StringVariant(ValueType val, unsigned int width)
    {
        std::ostringstream stream;
        stream << std::setprecision(width);
        stream << val;
        data.assign(stream.str());
    }

    template<typename ValueType>
    StringVariant& operator=(const ValueType val)
    {
        std::ostringstream stream;
        stream << val;
        data.assign(stream.str());

        return *this;
    }


    template<typename NumberType>
    NumberType toNumber() const
    {
        NumberType result = 0;
        std::istringstream stream(data);

        if(stream >> result)
            return result;
        else if(data == "yes" || data == "true")
            return 1;

        return 0;
    }

    bool toBool() const
    {
        if(data == "yes" || data == "true" || data == "t" || data == "1") return 1;
        return 0;
    }

    double toDouble() const
    {
        return toNumber<double>();
    }

    float toFloat() const
    {
        return toNumber<float>();
    }

    int toInt() const
    {
        return toNumber<int>();
    }

    std::string toString() const
    {
        return data;
    }

};

template<typename ValueType> std::string convertToString(ValueType val) { return StringVariant(val).toString(); }
template<typename ValueType> std::string convertToString(ValueType val, unsigned int width) { return StringVariant(val,width).toString(); }

}
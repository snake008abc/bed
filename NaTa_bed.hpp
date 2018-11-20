/**
 * @file BED.hpp
 * @brief a parser of BED format file
 *
 * @author JHH corp
 */

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/convert.hpp>
#include <boost/convert/lexical_cast.hpp>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace biovoltron::format                                          {
/// This is enumeration Bed_data for the function get_data<>()
enum Bed_data
{ 
        chr_name,
        chr_start,
        chr_end,
        name,
        score,
        strand,
        thi_start,
        thi_end,
        Rgb,
        blo_count,
        blo_size,
        blo_start
};


/**
 * @class Bed
 * @brief A class which stores Bed data in tuple type
 * 
 * @tparam TUPLETYPE: defaulted as tuple< string, uint32_t, 
 *  uint32_t>, for keeping track of the name of chromosome and its
 *  start and end positions
 * 
 * This class stores Bed entry and provides some functions for
 * construction, modification and entry access.
 */
template<class TUPLETYPE = std::tuple<std::string, int32_t, int32_t>>
class Bed
{
  private:
    ///A tuple storing the data
	TUPLETYPE data;

  public:
	///A default constructor with an empty tuple of data
	Bed():data()
	{}
	/**
	 * @brief A constructor with the current data
	 * 
	 * @param rhs: a lvalue reference of the tuple of data 
	 */
	Bed(TUPLETYPE& rhs):data(rhs)
	{}
	/**
	 * @brief A constructor by moving the known data
	 * 
	 * @param rhs: a rvalue reference of the tuple of data
	 */
	Bed(TUPLETYPE&& rhs):data(std::move(rhs))
	{}
	/**
	 * @brief A copy constructor
	 * 
	 * @param rhs: a lvalue reference of current Bed object
	 */
	Bed(const Bed& rhs):data(rhs.data)
	{}
	/**
	 * @brief A move constructor
	 * 
	 * @param rhs: a rvalue reference of the Bed object that user
	 *  want to copy
	 */
	Bed(Bed&& rhs):data(std::move(rhs.data))
	{}

	/**
	 * @brief operator = overload
	 * 
	 * @param rhs: the Bed object that user want to be assigned to
	 * @return Bed& 
	 */
	Bed& operator =(Bed& rhs)
	{
		data = rhs.data;
		return *this;
	}
	/**
	 * @brief operator = overload
	 * 
	 * @param rhs: the Bed object that user want to be assigned to
	 * @return Bed& 
	 */
	Bed& operator =(Bed&& rhs)
	{
		data = std::move(rhs.data);
		return *this;
	}
	
	///return Bed data in string format
	std::string to_string()
	{
		std::string rhs;
		get_string(rhs);
		return rhs;
	}

  private:
	
	/**
	 * @brief A function to transform all data saved in a tuple into
	 *  the string format. To be more specific, this function will
	 *  read and transform each element in tuple sequently.
	 * @tparam N: to record the index of the current element. Its
	 *  initial value is 0 and the function will end once N is larger
	 *  than or equal to the size of the original tuple.
	 * @param rhs: the string which is going to be written in the 
	 *  output string
	 */
	template<int N = 0>
	void get_string(std::string& rhs)
	{
		if constexpr (N >= std::tuple_size<TUPLETYPE>::value)
		{
			rhs.pop_back();
			return;
		}
		else
		{
			boost::cnv::lexical_cast cnv;
			boost::optional<std::string> v = boost::convert<std::string>(std::get<N>(data), cnv);
			rhs += (*v + '\t');
			get_string<N+1>(rhs);
		}
	}
  public:
	/**
	 * @brief << operator which gets entry from an ostream and it does
	 *  the same thing as function dump
	 * 
	 * @param rhs: An ostream where to dump entry
	 * @param i: A Bed object which the user want to dump
	 * @return std::ostream&: Is identical to parameter rhs
	 */
	friend std::ostream& operator <<(std::ostream& rhs, Bed& i)
	{
		rhs << i.to_string();
		return rhs;
	}
	/**
	 * @brief < opereator to show whether a Bed object is smaller 
	 *  than the other Bed object
	 * 
	 * @param b1: A Bed object the user want to compare
	 * @param b2: Another Bed object the user want to compare with b1
	 * @return true: If the data of b1 is smaller than the one of b2
	 * @return false: If the data of b1 is larger than or equals to 
	 *  the one of b2
	 */
	friend bool operator <(const Bed& b1, const Bed& b2)
	{
		return b1.data < b2.data;
	}

	/**
	 * @brief A function to dump a Bed object from a container to
	 *  an ostream
	 * 
	 * @param rhs: An ostream to dump Bed object
	 * @param c: A container of Bed
	 */
	static void dump(std::ostream& rhs, std::vector<Bed>& c)
	{
		std::string s;
		for(auto& i : c)
		{
			std::string t;
			t = i.to_string();
			s += (t+'\n');
		}
		rhs << s;
	}

	/**
	 * @brief >> operetor which gets entry from an istream and what
	 *  it does is basically same as function get_obj
	 * 
	 * @param rhs: An istream where to get entry 
	 * @param i: A Bed object to store entry
	 * @return std::istream&: Is identical to parameter is
	 */
	friend std::istream& operator >>(std::istream& rhs, Bed& i)
	{
		get_obj(rhs, i);
		return rhs;
	}

	/**
	 * @brief A function to load one chromosome from an istream.
	 *  Firstly, this function reads the content from an istream 
	 *  line-by-line and splits each line into a few of strings which
	 *  will be stored into a vector of string. Finally, the function
	 *  transforms this vector of string of each line into a tuple of
	 *  Bed object.
	 * 
	 * @param rhs: An istream contains the entry to data
	 * @param i: A Bed object to store all the data from istream
	 * @return std::istream&: Is identical to parameter rhs
	 */
	static std::istream& get_obj(std::istream& rhs, Bed& i)
	{
		std::string in, out;
		std::vector<std::string> v;
		if (std::getline(rhs, in))
		{
			boost::split (v, in, boost::is_any_of ("\t"));
			to_tuple(v, i);
		}
		return rhs;
	}

	/**
	 * @brief A funtion to transform a vector of string into a Bed
	 *  object. It will handle data separately to make sure they are 
	 *  in accordance with the format of each data in a Bed object.
	 * 
	 * @tparam N: to record the index of the current elements. Its
	 *  initial value is 0 and the function will end once N is larger
	 *  than or equal to the size of the original tuple.
	 * @param rhs: A vector of string in which the data in string
	 *  format is stored
	 * @param i: A Bed object which stores the data gained from a
	 *  vector of string
	 */
	template<int N = 0>
	static void to_tuple(std::vector<std::string>& rhs, Bed& i)
	{
		if constexpr (N >= std::tuple_size<TUPLETYPE>::value)
			return;
		else
		{
		/*
			using type = typename std::tuple_element<N, decltype(i.data)>::type;
			boost::cnv::lexical_cast cnv;
      boost::optional<type> pp = boost::convert<type>(rhs[N], cnv);
			std::get<N>(i.data) = *pp;
			to_tuple<N+1>(rhs, i);
		*/
			using type = typename std::tuple_element<N, decltype(i.data)>::type;
			if constexpr (std::is_same<int32_t, type>::value)
			{
				std::get<N>(i.data) = stoi(rhs[N]);
			}
			else if constexpr (std::is_same<std::string, type>::value)
			{
				std::get<N>(i.data) = rhs[N];
			}
			else if constexpr (std::is_same<char, type>::value)
			{
				std::get<N>(i.data) = rhs[N][0];
			}
			else
			{
				std::get<N>(i.data) = rhs[N];
			}
			to_tuple<N+1>(rhs, i);
		}
	}

	/**
	 * @brief An function to get the data object
	 * 
	 * @tparam N: the index of the desired data
	 * @return a: the data which the user want to get
	 */
	template<int N>
	auto get_data()
	{
		auto a = std::get<N>(data);
		return a;
	}

};

}



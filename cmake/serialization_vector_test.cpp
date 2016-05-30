
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <sstream>
#include <cstdio>

struct Int {
    Int(int _i=0) : i(_i) {};
    int i;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(i);
    }
};

int main()
{
    std::stringstream ss;
    {
        std::vector<Int> i;
        i.push_back(Int(10));
        boost::archive::xml_oarchive oa(ss);
        oa << BOOST_SERIALIZATION_NVP(i);
    }
    try
    {
        std::vector<Int> i;
        i.push_back(20);
        boost::archive::xml_iarchive ia(ss);
        ia >> BOOST_SERIALIZATION_NVP(i);
        if (i.size()>1)
            return 1;
    }
    catch (boost::archive::archive_exception &e)
    {
        return 1;
    }
    return 0;
}

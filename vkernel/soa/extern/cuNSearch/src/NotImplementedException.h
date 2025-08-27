#pragma once
#include <stdexcept>
#include <string>

class NotImplementedException : std::runtime_error
{
public:
	NotImplementedException()
		: std::runtime_error("Not implemented.")
	{

	}

	NotImplementedException(std::string str)
		: std::runtime_error("Not implemented: " + str)
	{
		std::cout << "Not implemented: " << str << std::endl;
	}
};
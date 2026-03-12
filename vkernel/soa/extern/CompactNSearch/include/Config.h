#pragma once

namespace CompactNSearch
{
#ifdef USE_DOUBLE
	using Real = double;
#else
	using Real = float;
#endif
}

#define INITIAL_NUMBER_OF_INDICES   50
#define INITIAL_NUMBER_OF_NEIGHBORS 50

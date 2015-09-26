#include "include_me.hpp"

template<class T1, class T2, class T3>
struct triple {
	T1 first;
	T2 second;
	T3 third;

	triple() :
		first(),
		second(),
		third() { }

	triple(T1 new_first, T2 new_second, T3 new_third) :
		first(new_first),
		second(new_second),
		third(new_third) { }
};

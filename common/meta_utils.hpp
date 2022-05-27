#ifndef __META_UTILS_HPP__
#define __META_UTILS_HPP__

#include <type_traits>

/******************************************
 *   Meta Utils (v0.1)                    *
 *       Some useful meta programming     *
 *           utils written in c++11.      *
 *                                        *
 *       Author: zexlus1126               *
 *       Date: 2018/05/09                 *
 ******************************************/

namespace zex
{

// ===== is_any =====
//
// tempate<class T, class V1, class V2, ...>
//
// Trait class that identifies whether T is the same type 
// as any type in the (V1, V2 ...) list
//
// Usage:
//   int main()
//   {
//      std::cout << std::boolalpha;
//      std::cout << zex::is_any<int, char, 
//                                    short, 
//                                    double, 
//                                    std::vector<int>, 
//                                    std::string, 
//                                    int
//                                    >::value << std::endl;
//      return 0;
//   }
// 
// Output: 
//    true
// Description:
//    Type int is same as the 6th element contains in the type list.

template<typename T, typename ... R>
struct is_any : std::false_type {};

template<typename T, typename F>
struct is_any<T, F> : std::is_same<T, F> {};

template<typename T, typename F, typename ...R>
struct is_any<T, F, R...> :
        std::integral_constant<bool, std::is_same<T, F>::value
                                    || is_any<T, R...>::value>{};


// ===== type_index =====
// 
// template<class T, class V1, class V2, ...>
//
// Trait class that identifies whether T is the same type
// as any type in the (V1, V2 ...) list, and then return
// the index number of the first matched element.
//
// Usage:
//   int main()
//   {
//      std::cout << zex::type_index<int, char,    //0
//                                        short,   //1
//                                        double,  //2
//                                        int      //3
//                                        >::value << std::endl;
//      return 0;
//   }
//
// Output:
//    3

template<int _case, typename T, typename ...R>
struct _type_index : std::false_type {};

template<int _case, typename T, typename F>
struct _type_index<_case, T, F> : std::conditional<std::is_same<T, F>::value,
                                        std::integral_constant<int, _case>,
                                        std::integral_constant<int, -1>>::type{};

template<int _case, typename T, typename F, typename ...R>
struct _type_index<_case, T, F, R...> : std::conditional<std::is_same<T, F>::value,
                                        _type_index<_case, T, F>, _type_index<_case+1, T, R...>>::type{};

template<typename T, typename ...R>
struct type_index : _type_index<0, T, R...>{};



// ===== type_case =====
//
// template<int case, class V1, class V2, ...>
//
// Obtions the specified class type, depending on the first 
// template parameter 'case' spesifing the index of the 
// element in the type list (V1, V2, ...).
//
// Usage:
//
//   template<int T>
//   using tcase = zex::type_case<T, int,        //0
//                                   unsigned,   //1
//                                   double,     //2
//                                   >;
//   int main()
//   {
//
//      std::cout << std::boolalpha;
//      std::cout << std::is_same<tcase<2>::type, double>::value << std::endl;
//      return 0;
//   }
//
// Output:
//   true

template<int _case, typename ...R>
struct _type_case : std::false_type{};

template<int _case, typename T>
struct _type_case<_case, T> {  typedef T type;  };

template<int _case, typename T, typename ...R>
struct _type_case<_case, T, R...> : std::conditional<_case==0, _type_case<_case, T>,
                                                               _type_case<_case-1, R...>>::type{};

template<int _case, typename ...T>
struct type_case : _type_case<_case, T...> {};


// ===== tv_pair =====
//
// template<typename T, T V>
//
// This class wrap the type and the constant value together
//
// Usage:
//
//   int main()
//   {
//      std::cout << "int: " << zex::tv_pair<int, 10>::value << std::endl
//                << "char: " << zex::tv_pair<char, 'a'>::value << std::endl;
//      return 0;
//   }
//
// Output:
//   10
//   a

template<typename T, T V>
struct tv_pair {
    typedef T type;
    static constexpr T value=V;
};


// ===== pick/options =====

template<int _case>
struct pick{
    constexpr static bool _options() { return false; }

    template<typename T, typename ...Tail>
    constexpr static auto _options(const T& head, Tail&&...tail)
        -> decltype(pick<_case-1>::_options(std::forward<Tail>(tail)...))
    {
        return pick<_case-1>::_options(tail...);
    }

    template<typename T, typename ...Tail>
    constexpr static auto options(const T& head, Tail&&...tail)
        -> decltype(pick<_case-1>::_options(tail...))
    {
        return pick<_case-1>::_options(tail...);
    }
};

template<>
struct pick<0>{
    template<typename T, typename ...Tail>
    constexpr static auto _options(const T& head, Tail&&...tail)
        -> decltype(head)
    {
        return head;
    }

    template<typename T, typename ...Tail>
    constexpr static auto options(const T& head, Tail&&...tail)
        -> decltype(head)
    {
        return head;
    }
};

// namespace zex {end}
}


#endif

// Copyright (C) 2019 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_STL_INTERFACES_FWD_HPP
#define BOOST_STL_INTERFACES_FWD_HPP

#include <iterator>

#ifndef BOOST_STL_INTERFACES_DOXYGEN

// It is not known whether intel has fixed this issue, but it was
// observed to be present in icc 19.1.2.20200623. Update the version
// check when it is fixed.
#if defined(_MSC_VER) || defined(__GNUC__) && __GNUC__ < 8 || defined(__INTEL_COMPILER) && __INTEL_COMPILER < 1920
#define BOOST_STL_INTERFACES_NO_HIDDEN_FRIEND_CONSTEXPR
#define BOOST_STL_INTERFACES_HIDDEN_FRIEND_CONSTEXPR
#else
#define BOOST_STL_INTERFACES_HIDDEN_FRIEND_CONSTEXPR constexpr
#endif

#if defined(__GNUC__) && __GNUC__ < 9
#define BOOST_STL_INTERFACES_CONCEPT concept bool
#else
#define BOOST_STL_INTERFACES_CONCEPT concept
#endif

#endif


namespace gmx { namespace boost { namespace stl_interfaces {
    inline namespace v1 {

        /** An enumeration used to indicate whether the underlying data have a
            contiguous or discontiguous layout when instantiating
            `view_interface` and `sequence_container_interface`. */
        enum class element_layout : bool {
            discontiguous = false,
            contiguous = true
        };

        namespace v1_dtl {
            template<typename... T>
            using void_t = void;

            template<typename Iter>
            using iter_difference_t =
                typename std::iterator_traits<Iter>::difference_type;

            template<typename Range, typename = void>
            struct iterator;
            template<typename Range>
            struct iterator<
                Range,
                void_t<decltype(std::declval<Range &>().begin())>>
            {
                using type = decltype(std::declval<Range &>().begin());
            };
            template<typename Range>
            using iterator_t = typename iterator<Range>::type;

            template<typename Range, typename = void>
            struct sentinel;
            template<typename Range>
            struct sentinel<
                Range,
                void_t<decltype(std::declval<Range &>().end())>>
            {
                using type = decltype(std::declval<Range &>().end());
            };
            template<typename Range>
            using sentinel_t = typename sentinel<Range>::type;

            template<typename Range>
            using range_difference_t = iter_difference_t<iterator_t<Range>>;

            template<typename Range>
            using common_range =
                std::is_same<iterator_t<Range>, sentinel_t<Range>>;

            template<typename Range, typename = void>
            struct decrementable_sentinel : std::false_type
            {
            };
            template<typename Range>
            struct decrementable_sentinel<
                Range,
                void_t<decltype(--std::declval<sentinel_t<Range> &>())>>
                : std::true_type
            {
            };
        }

    }
}}}

#endif

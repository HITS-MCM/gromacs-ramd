/*
 * vectypes_io.h
 *
 *  Created on: Jun 25, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include <iostream>

#include "vectypes.h"

/// Print BasicVector
template <typename T>
std::ostream& operator << (std::ostream& os, gmx::BasicVector<T> const& v)
{
    return os << "[" << v[0] << " ," << v[1] << " ," << v[2] << "]";
}

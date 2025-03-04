/*
 * Copyright 2022 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(multi_UAP.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(8a14b49942d645cf1ef9dd17587c7da8)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <bluetooth/multi_UAP.h>
// pydoc.h is automatically generated in the build directory
#include <multi_UAP_pydoc.h>

void bind_multi_UAP(py::module& m)
{

    using multi_UAP    = ::gr::bluetooth::multi_UAP;


    py::class_<multi_UAP, gr::bluetooth::multi_block,
        std::shared_ptr<multi_UAP>>(m, "multi_UAP", D(multi_UAP))

        .def(py::init(&multi_UAP::make),
           py::arg("sample_rate"),
           py::arg("center_freq"),
           py::arg("squelch_threshold"),
           py::arg("LAP"),
           D(multi_UAP,make)
        )
        



        ;




}









/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __QUAD_QUADIRON_H__
#define __QUAD_QUADIRON_H__

/** @namespace quadiron
 *
 * The root namespace of the QuadIron library.
 */

#include "build_info.h"
#include "fec_base.h"
#include "fec_rs_fnt.h"
#include "fec_rs_gf2n.h"
#include "fec_rs_gf2n_fft.h"
#include "fec_rs_gf2n_fft_add.h"
#include "fec_rs_gfp_fft.h"
#include "fec_rs_nf4.h"

/** Return the version string of QuadIron.
 *
 * The version string has the form MAJOR.MINOR.PATCH-REVISION, where '-REVISION'
 * is optional (only present for development version).
 *
 * @return the version string.
 */
static inline const char* get_version()
{
    return VERSION;
}

#endif

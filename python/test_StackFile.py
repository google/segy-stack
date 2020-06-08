#!/usr/bin/python
#
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import pytest
import hashlib

from segystack import StackFile, SegyFile
from segystack.test import create_test_segy


def hash_data(data):
    h = hashlib.md5()
    h.update(data.tostring())
    return h.hexdigest()


class TestStackFile:

    def test_basic(self):
        with pytest.raises(TypeError):
            StackFile()

        opts = StackFile.SegyOptions()
        with pytest.raises(TypeError):
            opts.set_utm_zone(35, "A")
        with pytest.raises(TypeError):
            opts.set_utm_zone(120, "H")

    @pytest.mark.parametrize(
        "num_il,num_xl,num_ds",
        [
            pytest.param(
                101, 201, 601
            ),
            pytest.param(
                99, 199, 299
            ),
            pytest.param(
                1, 201, 601
            ),
            pytest.param(
                201, 1, 601
            ),
            pytest.param(
                301, 51, 1
            ),
        ],
    )
    def test_init(self, num_il, num_xl, num_ds):
        opts = StackFile.SegyOptions()
        opts.set_utm_zone(34, "H")

        sgy_file = "/tmp/test_segystack_" + \
            str(num_il) + "_" + str(num_xl) + "_" + str(num_ds) + ".sgy"
        samp_int = 1 + int(2000 / num_ds)
        il_incr = 1 + int(num_ds / num_il)
        xl_incr = 1 + int(num_ds / num_xl)
        il_sp = 20.0 * il_incr
        xl_sp = 10.0 * xl_incr
        create_test_segy(sgy_file, num_ds, samp_int, num_il,
                         il_incr, num_xl, xl_incr, 0.0, 0.0, il_sp, xl_sp, opts)

        sgy = SegyFile(sgy_file)
        sgy.open("r")
        offsets = sgy.guess_trace_header_offsets()
        opts.set_trace_header_offsets(offsets)

        outfile = "/tmp/out.stack"
        sf = StackFile(outfile, sgy, opts)
        sgy.close()

        self.check_stackfile(sf, opts, outfile, num_il,
                             num_xl, num_ds, samp_int, il_incr, xl_incr, il_sp, xl_sp)

        sf_in = StackFile(outfile)
        self.check_stackfile(sf_in, opts, outfile, num_il,
                             num_xl, num_ds, samp_int, il_incr, xl_incr, il_sp, xl_sp)

        os.remove(sgy_file)
        os.remove(outfile)
        os.remove(outfile + "_data")

    def check_stackfile(self, sf, opts, outfile, num_il, num_xl, num_ds,
                        samp_int, il_incr, xl_incr, il_sp, xl_sp):

        if (sf.grid().num_inlines > 1):
            assert(sf.grid().inline_increment == il_incr)
            assert(sf.grid().inline_spacing == il_sp)

        il_nums = set(sf.grid().inline_numbers)
        for il in range(sf.grid().inline_min, sf.grid().inline_max+1, sf.grid().inline_increment):
            assert(il in il_nums)
        with pytest.raises(TypeError):
            sf.grid().inline_numbers = (10, 20)

        if (sf.grid().num_crosslines > 1):
            assert(sf.grid().crossline_increment == xl_incr)
            assert(sf.grid().crossline_spacing == xl_sp)

        xl_nums = set(sf.grid().crossline_numbers)
        for xl in range(sf.grid().crossline_min, sf.grid().crossline_max+1, sf.grid().crossline_increment):
            assert(xl in xl_nums)
        with pytest.raises(TypeError):
            sf.grid().crossline_numbers = (10, 20)

        utm_zone = sf.grid().utm_zone()
        assert(utm_zone.letter == opts.utm_zone().letter)
        assert(utm_zone.number == opts.utm_zone().number)
        assert(sf.grid().num_inlines == num_il)
        assert(sf.grid().num_crosslines == num_xl)
        assert(sf.grid().num_samples == num_ds)
        assert(abs(sf.grid().sampling_interval - samp_int/1000.0) < 1e-4)

        bbox = sf.grid().bounding_box()
        assert(bbox.c1.x == 0.0)
        assert(bbox.c1.y == 0.0)
        assert(bbox.c1.inline_num == 0)
        assert(bbox.c1.crossline_num == 0)
        assert(bbox.c1.latitude != 0.0)
        assert(bbox.c1.longitude != 0.0)
        if (sf.grid().num_inlines > 1):
            assert(bbox.c3.y != 0.0)
            assert(bbox.c4.y != 0.0)
            assert(bbox.c3.inline_num == (num_il-1)*il_incr)
            assert(bbox.c4.inline_num == (num_il-1)*il_incr)

        if (sf.grid().num_crosslines > 1):
            assert(bbox.c2.x != 0.0)
            assert(bbox.c4.x != 0.0)
            assert(bbox.c2.crossline_num == (num_xl-1)*xl_incr)
            assert(bbox.c4.crossline_num == (num_xl-1)*xl_incr)

        il_data = sf.read_inline(
            sf.grid().inline_numbers[int(sf.grid().num_inlines/2)], 0.0)
        assert(il_data.shape == (num_xl, num_ds))

        sf.set_crossline_access_opt(True)
        xl_data1 = sf.read_crossline(
            sf.grid().crossline_numbers[int(sf.grid().num_crosslines/2)], 0.0)
        assert(xl_data1.shape == (num_il, num_ds))
        sf.set_crossline_access_opt(False)

        with pytest.raises(IOError):
            open(outfile + "_data_xline", "r")
        xl_data2 = sf.read_crossline(
            sf.grid().crossline_numbers[int(sf.grid().num_crosslines/2)], 0.0)
        assert(xl_data2.shape == (num_il, num_ds))

        assert(hash_data(xl_data1) == hash_data(xl_data2))

        sf.set_depth_slice_access_opt(True)
        assert(open(outfile + "_data_depth", "r"))

        ds_data1 = sf.read_depth_slice(int(num_ds / 2), 0.0)
        assert(ds_data1.shape == (num_il, num_xl))

        sf.set_depth_slice_access_opt(False)
        with pytest.raises(IOError):
            open(outfile + "_data_depth", "r")

        ds_data2 = sf.read_depth_slice(int(num_ds / 2), 0.0)
        assert(ds_data2.shape == (num_il, num_xl))

        assert(hash_data(ds_data1) == hash_data(ds_data2))


if __name__ == "__main__":
    t = TestStackFile()
    t.test_basic()
    t.test_init()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate Observation Blocks for the 'Awakening the Giants' program on the
Wide Field Imager at the 2.2 m MPG telescope at La Silla Observatory.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import astropy.units as u
import astropy.coordinates as coord

import ob_templates

# Produce a set of calibration OBs, including a focus OB

# 10 Sky flats
# 10 Dome flats
# Biases?
# Focus sequence

clobber = True
ins_filter_name = "FILTERNAME"
exposure_time = 200




galactic_field_positions = [
    (-4, 0),
    (-4, -2),
    (-4, +2),
    (0, 0),
    (0, -2),
    (0, +2),
    (+4, 0),
    (+4, -2),
    (+4, +2)
]

repr_gal_deg = lambda x: "{0}{1:0>2d}".format("mp"[x > 0], abs(int(x * 10)))
generate_ob_name = lambda l, b: "atg.sci.{l_str}{b_str}".format(
    l_str=repr_gal_deg(l), b_str=repr_gal_deg(b))

N = len(galactic_field_positions)
for i, (l, b) in enumerate(galactic_field_positions):

    coordinate = coord.Galactic(l, b, unit=u.degree).transform_to(coord.ICRS)

    # Create the observation block.
    observation_block = ob_templates.ScienceObservationBlock(
        ob_name=generate_ob_name(l, b),
        ins_filter_name=ins_filter_name,
        target_ra=coordinate.ra,
        target_dec=coordinate.dec,
        exposure_time=exposure_time, 
    )

    observation_block.write(
        "{}.obx".format(observation_block.ob_name), clobber=clobber)

    logging.info(
        "Generated OB {0}/{1}: {2}".format(i, N, observation_block.ob_name))

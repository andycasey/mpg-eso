#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Templates of Observation Blocks for the Wide Field Imager.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import os
import yaml


class ObservationBlock(object):

    _default_path = os.path.join(os.path.dirname(__file__), "defaults.yaml")
    _template_path = None

    def __init__(self, ob_name, ins_filter_name, user_comments=None, **kwargs):
        """
        An abstract observation block object.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        # Load in the defaults.
        with open(self._default_path, "r") as fp:
            defaults = yaml.load(fp)
        self.__dict__.update(**defaults)

        self.__dict__.update({
            "ob_name": ob_name,
            "observation_description_name": ob_name,
            "user_comments": user_comments or "",
            "ins_filter_name": ins_filter_name,
        })

        return None


    def write(self, filename, clobber=False):
        """
        Write the observing block to disk.

        :param filename:
            The path to write the filename to.

        :param clobber: [optional]
            Overwrite the filename if it already exists.
        """

        assert self._template_path is not None, \
            "class _template_path cannot be None"

        if os.path.exists(filename) and not clobber:
            raise IOError("filename {} already exists".format(filename))

        if not filename.lower().endswith(".obx"):
            logging.warn("Filename '{}' does not have an OBX extension".format(
                filename))

        # Load and format the contents.
        with open(self._template_path, "r") as fp:
            template = fp.read()

        contents = template.format(**self.__dict__)

        with open(filename, "w") as fp:
            fp.write(contents)

        return True



class FlatObservationBlock(ObservationBlock):

    _template_path = \
        os.path.join(os.path.dirname(__file__), "flat.template.obx")

    def __init__(self, ob_name, ins_filter_name, user_comments=None, **kwargs):
        """
        An abstract flat field observation block.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        super(FlatObservationBlock, self).__init__(
            ob_name, ins_filter_name, user_comments, **kwargs)
        return None


class SkyFlatObservationBlock(FlatObservationBlock):

    def __init__(self, ob_name, ins_filter_name, user_comments=None, **kwargs):
        """
        A sky flat observation block.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        super(SkyFlatObservationBlock, self).__init__(
            ob_name, ins_filter_name, user_comments, **kwargs)

        self.__dict__.update({
            "template_name": "WFI_img_cal_SkyFlat"
        })
        return None


class DomeFlatObservationBlock(FlatObservationBlock):

    def __init__(self, ob_name, ins_filter_name, user_comments=None, **kwargs):
        """
        A dome flat observation block.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        super(DomeFlatObservationBlock, self).__init__(
            ob_name, ins_filter_name, user_comments, **kwargs)

        self.__dict__.update({
            "template_name": "WFI_img_cal_DomeFlat"
        })
        return None


class FocusSequenceObservationBlock(ObservationBlock):

    _template_path = \
        os.path.join(os.path.dirname(__file__), "focus.template.obx")


    def __init__(self, ob_name, ins_filter_name, user_comments=None, **kwargs):
        """
        A focus sequence observation block.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        super(FocusSequenceObservationBlock, self).__init__(
            ob_name, ins_filter_name, user_comments, **kwargs)

        self.__dict__.update({
            "template_name": "WFI_cal_FocusSeq",
            "det_win1_uit1": 15,
            "det_win1_offset": 50,
            "seq_nexpo": 9,
        })
        return True


class ScienceObservationBlock(ObservationBlock):

    _template_path = \
        os.path.join(os.path.dirname(__file__), "science.template.obx")

    def __init__(self, ob_name, ins_filter_name, exposure_time, 
        user_comments=None, **kwargs):
        """
        A science object observation block.

        :param ob_name:
            The name to assign of the Observation Block.

        :param ins_filter_name:
            The name of the instrument filter to use.

        :param user_comments: [optional]
            Any user comments to specify for the Observation Block.
        """

        super(ScienceObservationBlock, self).__init__(
            ob_name, ins_filter_name, user_comments, **kwargs)

        self.__dict__.update({
            "det_win1_uit1": exposure_time,
            "seq_nexpo": 1
        })

        return True

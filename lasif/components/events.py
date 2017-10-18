#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import copy
import glob
import os
import warnings

import obspy
import pyasdf
from lasif import LASIFNotFoundError, LASIFWarning
from obspy.geodetics import FlinnEngdahl

from .component import Component


class EventsComponent(Component):
    """
    Component managing a folder of QuakeML files.

    Each file must adhere to the scheme ``*.xml``.

    :param folder: Folder with QuakeML files.
    :param communicator: The communicator instance.
    :param component_name: The name of this component for the
        communicator.
    """

    def __init__(self, folder, communicator, component_name):
        super(EventsComponent, self).__init__(communicator, component_name)
        self.__event_info_cache = {}
        self.folder = folder

        self.index_values = [
            ("filename"),
            ("event_name"),
            ("latitude"),
            ("longitude"),
            ("depth_in_km"),
            ("magnitude"),
            ("region")]

        self.all_events = {}
        self.fill_all_events()

    def fill_all_events(self):
        files = glob.glob(os.path.join(self.folder, '*.h5'))
        for file in files:
            event_name = os.path.splitext(os.path.basename(file))[0]
            self.all_events[event_name] = file

    def update_cache(self):
        files = glob.glob(os.path.join(self.folder, '*.h5'))
        for file in files:
            event_name = os.path.splitext(os.path.basename(file))[0]
            self.get(event_name)

    @staticmethod
    def _extract_index_values_quakeml(filename):
        """
        Reads QuakeML files and extracts some keys per channel. Only one
        event per file is allowed.
        """
        ds = pyasdf.ASDFDataSet(filename, mode='r')

        event_name = os.path.splitext(os.path.basename(filename))[0]
        org = ds.waveforms[event_name].coordinates

        return [
            str(filename),
            str(event_name),
            float(org["latitude"]),
            float(org["longitude"]),
            float(1.0/1000.0),
            float(1.0e10),
            str(FlinnEngdahl().get_region(org["longitude"], org["latitude"]))
        ]

    def list(self):
        """
        List of all events.
        """
        self.update_cache()
        return sorted(self.__event_info_cache.keys())

    def count(self):
        """
        Get the number of events managed by this component.
        """
        return len(self.all_events)

    def has_event(self, event_name):
        """
        Test for existence of an event.

        :type event_name: str
        :param event_name: The name of the event.
        """
        # Make sure it also works with existing event dictionaries. This
        # has the potential to simplify lots of code.
        try:
            event_name = event_name["event_name"]
        except (KeyError, TypeError):
            pass
        return event_name in self.all_events

    def get_all_events(self):
        """
        Returns a dictionary with the key being the event names and the
        values the information about each event, as would be returned by the
        :meth:`~lasif.components.events.EventsComponent.get` method.
        """
        # make sure cache is filled
        self.update_cache()
        return copy.deepcopy(self.__event_info_cache)

    def get(self, event_name):
        """
        Get information about one event.

        This function uses multiple cache layers and is thus very cheap to
        call.

        :type event_name: str
        :param event_name: The name of the event.
        :rtype: dict
        """
        try:
            event_name = event_name["event_name"]
        except (KeyError, TypeError):
            pass

        if event_name not in self.all_events:
            raise LASIFNotFoundError("Event '%s' not known to LASIF." %
                                     event_name)

        if event_name not in self.__event_info_cache:
            values = dict(zip(self.index_values,
                              self._extract_index_values_quakeml(
                                  self.all_events[event_name])))
            self.__event_info_cache[event_name] = values
        return self.__event_info_cache[event_name]

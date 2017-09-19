#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import copy
import warnings

import os
import glob

import obspy

from .component import Component
from lasif import LASIFNotFoundError, LASIFWarning
from obspy.geodetics import FlinnEngdahl
import pyasdf


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
            ("origin_time"),
            ("m_rr"),
            ("m_pp"),
            ("m_tt"),
            ("m_rp"),
            ("m_rt"),
            ("m_tp"),
            ("magnitude"),
            ("magnitude_type"),
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
        for filename in files:
            event_name = os.path.splitext(os.path.basename(filename))[0]
            self.get(event_name)

    @staticmethod
    def _extract_index_values_quakeml(filename):
        """
        Reads QuakeML files and extracts some keys per channel. Only one
        event per file is allows.
        """
        ds = pyasdf.ASDFDataSet(filename, mode='r')
        event = ds.events[0]

        # Extract information.
        mag = event.preferred_magnitude() or event.magnitudes[0]
        org = event.preferred_origin() or event.origins[0]
        if org.depth is None:
            warnings.warn("Origin contains no depth. Will be assumed to be 0",
                          LASIFWarning)
            org.depth = 0.0
        if mag.magnitude_type is None:
            warnings.warn("Magnitude has no specified type. Will be assumed "
                          "to be Mw", LASIFWarning)
            mag.magnitude_type = "Mw"

        # Get the moment tensor.
        fm = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
        mt = fm.moment_tensor.tensor

        event_name = os.path.splitext(os.path.basename(filename))[0]

        return [
            str(filename),
            str(event_name),
            float(org.latitude),
            float(org.longitude),
            float(org.depth / 1000.0),
            float(org.time.timestamp),
            float(mt.m_rr),
            float(mt.m_pp),
            float(mt.m_tt),
            float(mt.m_rp),
            float(mt.m_rt),
            float(mt.m_tp),
            float(mag.mag),
            str(mag.magnitude_type),
            str(FlinnEngdahl().get_region(org.longitude, org.latitude))
        ]

    def list(self):
        """
        List of all events.

        >>> comm = getfixture('events_comm')
        >>> comm.events.list() #  doctest: +NORMALIZE_WHITESPACE
        ['GCMT_event_ICELAND_Mag_5.5_2014-10-7-10',
        'GCMT_event_IRAN-IRAQ_BORDER_REGION_Mag_5.8_2014-10-15-13']
        """
        self.update_cache()
        return sorted(self.__event_info_cache.keys())

    def count(self):
        """
        Get the number of events managed by this component.

        >>> comm = getfixture('events_comm')
        >>> comm.events.count()
        2
        """
        return len(self.all_events)

    def has_event(self, event_name):
        """
        Test for existence of an event.

        :type event_name: str
        :param event_name: The name of the event.

        >>> comm = getfixture('events_comm')
        >>> comm.events.has_event('GCMT_event_ICELAND_Mag_5.5_2014-10-7-10')
        True
        >>> comm.events.has_event('random')
        False
        """
        # Make sure  it also works with existing event dictionaries. This
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

        >>> comm = getfixture('events_comm')
        >>> comm.events.get_all_events() \
        # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        {'GCMT_event_ICELAND_Mag_5.5_2014-10-7-10': {...},
         'GCMT_event_IRAN-IRAQ_BORDER_REGION_Mag_5.8_2014-10-15-13': {...}}
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

        >>> comm = getfixture('events_comm')
        >>> comm.events.get('GCMT_event_ICELAND_Mag_5.5_2014-10-7-10') \
        # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        {'filename': '...',
        'event_name': 'GCMT_event_ICELAND_Mag_5.5_2014-10-7-10',
        'latitude': 64.62,
        'longitude': -17.26, 'depth_in_km': 12.0, 'origin_time':
        UTCDateTime(2014, 10, 7, 10, 22, 34, 100000), 'm_rr': -2.97e+17,
         'm_pp': 1.23e+17,
        'm_tt': 1.74e+17, 'm_rp': -7.74e+16, 'm_rt': 3.87e+16, 'm_tp':
        -1900000000000000.0, 'magnitude': 5.53, 'magnitude_type': 'Mwc',
         'region': 'ICELAND'}

        The moment tensor components are in ``Nm``. The dictionary will
        contain the following keys:

        >>> sorted(comm.events.get(
        ...     'GCMT_event_ICELAND_Mag_5.5_2014-10-7-10').keys()) \
        # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
         ['depth_in_km', 'event_name', 'filename', 'latitude', 'longitude',
          'm_pp', 'm_rp', 'm_rr', 'm_rt', 'm_tp', 'm_tt', 'magnitude',
          'magnitude_type', 'origin_time', 'region']

        It also works with an existing event dictionary. This eases calling
        the function under certain circumstances.

        >>> ev = comm.events.get('GCMT_event_ICELAND_Mag_5.5_2014-10-7-10')
        >>> ev == comm.events.get(ev)
        True
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
            values["origin_time"] = obspy.UTCDateTime(values["origin_time"])
            self.__event_info_cache[event_name] = values
        return self.__event_info_cache[event_name]

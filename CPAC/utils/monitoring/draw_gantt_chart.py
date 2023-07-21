# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Resolves bugs preventing the original from generating the chart
#     * Handles when chart-drawing is called but no nodes were run (all cached)

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#     Copyright (c) 2015-2019, Nipype developers

#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at

#         http://www.apache.org/licenses/LICENSE-2.0

#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

#     Prior to release 0.12, Nipype was licensed under a BSD license.

# Modifications Copyright (C) 2021-2023 C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
'''Module to draw an html gantt chart from logfile produced by
``CPAC.utils.monitoring.log_nodes_cb()``.

See https://nipype.readthedocs.io/en/latest/api/generated/nipype.utils.draw_gantt_chart.html
'''  # noqa: E501
import random

from collections import OrderedDict
from datetime import datetime
from warnings import warn

from nipype.utils.draw_gantt_chart import draw_lines, draw_resource_bar, \
                                          log_to_dict


def create_event_dict(start_time, nodes_list):
    """
    Function to generate a dictionary of event (start/finish) nodes
    from the nodes list

    Parameters
    ----------
    start_time : datetime.datetime
        a datetime object of the pipeline start time
    nodes_list : list
        a list of the node dictionaries that were run in the pipeline

    Returns
    -------
    events : dictionary
        a dictionary where the key is the timedelta from the start of
        the pipeline execution to the value node it accompanies
    """

    # Import packages
    import copy

    events = {}
    for node in nodes_list:
        # Format node fields
        estimated_threads = node.get("num_threads", 1)
        estimated_memory_gb = node.get("estimated_memory_gb", 1.0)
        runtime_threads = node.get("runtime_threads", 0)
        runtime_memory_gb = node.get("runtime_memory_gb", 0.0)

        # Init and format event-based nodes
        node["estimated_threads"] = estimated_threads
        node["estimated_memory_gb"] = estimated_memory_gb
        node["runtime_threads"] = runtime_threads
        node["runtime_memory_gb"] = runtime_memory_gb
        start_node = node
        finish_node = copy.deepcopy(node)
        start_node["event"] = "start"
        finish_node["event"] = "finish"

        # Get dictionary key
        start_delta = (node["start"] - start_time).total_seconds()
        finish_delta = (node["finish"] - start_time).total_seconds()

        # Populate dictionary
        if events.get(start_delta):
            err_msg = "Event logged twice or events started at exact same " \
                      "time!"
            warn(str(KeyError(err_msg)), category=Warning)
        events[start_delta] = start_node
        events[finish_delta] = finish_node

    # Return events dictionary
    return events


def calculate_resource_timeseries(events, resource):
    """
    Given as event dictionary, calculate the resources used
    as a timeseries

    Parameters
    ----------
    events : dictionary
        a dictionary of event-based node dictionaries of the workflow
        execution statistics
    resource : string
        the resource of interest to return the time-series of;
        e.g. 'runtime_memory_gb', 'estimated_threads', etc

    Returns
    -------
    time_series : pandas Series
        a pandas Series object that contains timestamps as the indices
        and the resource amount as values
    """
    # Import packages
    import pandas as pd

    # Init variables
    res = OrderedDict()
    all_res = 0.0

    # Iterate through the events
    for _, event in sorted(events.items()):
        if event["event"] == "start":
            if resource in event:
                try:
                    all_res += float(event[resource])
                except ValueError:
                    next
            current_time = event["start"]
        elif event["event"] == "finish":
            if resource in event:
                try:
                    all_res -= float(event[resource])
                except ValueError:
                    next
            current_time = event["finish"]
        res[current_time] = all_res

    # Formulate the pandas timeseries
    time_series = pd.Series(data=list(res.values()), index=list(res.keys()))
    # Downsample where there is only value-diff
    ts_diff = time_series.diff()
    time_series = time_series[ts_diff != 0]

    # Return the new time series
    return time_series


def draw_nodes(start, nodes_list, cores, minute_scale, space_between_minutes,
               colors):
    """
    Function to return the html-string of the node drawings for the
    gantt chart

    Parameters
    ----------
    start : datetime.datetime obj
        start time for first node
    nodes_list : list
        a list of the node dictionaries
    cores : integer
        the number of cores given to the workflow via the 'n_procs'
        plugin arg
    total_duration : float
        total duration of the workflow execution (in seconds)
    minute_scale : integer
        the scale, in minutes, at which to plot line markers for the
        gantt chart; for example, minute_scale=10 means there are lines
        drawn at every 10 minute interval from start to finish
    space_between_minutes : integer
        scale factor in pixel spacing between minute line markers
    colors : list
        a list of colors to choose from when coloring the nodes in the
        gantt chart

    Returns
    -------
    result : string
        the html-formatted string for producing the minutes-based
        time line markers
    """
    # Init variables
    result = ""
    scale = space_between_minutes / minute_scale
    space_between_minutes = space_between_minutes / scale
    end_times = [
        datetime(
            start.year, start.month, start.day, start.hour, start.minute,
            start.second
        )
        for core in range(cores)
    ]

    # For each node in the pipeline
    for node in nodes_list:
        # Get start and finish times
        node_start = node["start"]
        node_finish = node["finish"]
        # Calculate an offset and scale duration
        offset = (
            (node_start - start).total_seconds() / 60
        ) * scale * space_between_minutes + 220
        # Scale duration
        scale_duration = (node["duration"] / 60) * scale \
            * space_between_minutes
        if scale_duration < 5:
            scale_duration = 5
        scale_duration -= 2
        # Left
        left = 60
        for core in range(len(end_times)):
            if end_times[core] < node_start:
                left += core * 30
                end_times[core] = datetime(
                    node_finish.year,
                    node_finish.month,
                    node_finish.day,
                    node_finish.hour,
                    node_finish.minute,
                    node_finish.second,
                )
                break

        # Get color for node object
        color = random.choice(colors)
        if "error" in node:
            color = "red"

        # Setup dictionary for node html string insertion
        node_dict = {
            "left": left,
            "offset": offset,
            "scale_duration": scale_duration,
            "color": color,
            "node_name": node.get("name", node.get("id", "")),
            "node_dur": node["duration"] / 60.0,
            "node_start": node_start.strftime("%Y-%m-%d %H:%M:%S"),
            "node_finish": node_finish.strftime("%Y-%m-%d %H:%M:%S"),
        }
        # Create new node string
        new_node = (
            "<div class='node' style='left:%(left)spx;top:%(offset)spx;"
            "height:%(scale_duration)spx;background-color:%(color)s;'"
            "title='%(node_name)s\nduration:%(node_dur)s\n"
            "start:%(node_start)s\nend:%(node_finish)s'></div>" % node_dict
        )

        # Append to output result
        result += new_node

    # Return html string for nodes
    return result


def generate_gantt_chart(
    logfile,
    cores,
    minute_scale=10,
    space_between_minutes=50,
    colors=["#7070FF", "#4E4EB2", "#2D2D66", "#9B9BFF"],
):
    """
    Generates a gantt chart in html showing the workflow execution based on a callback log file.
    This script was intended to be used with the MultiprocPlugin.
    The following code shows how to set up the workflow in order to generate the log file:

    Parameters
    ----------
    logfile : string
        filepath to the callback log file to plot the gantt chart of
    cores : integer
        the number of cores given to the workflow via the 'n_procs'
        plugin arg
    minute_scale : integer (optional); default=10
        the scale, in minutes, at which to plot line markers for the
        gantt chart; for example, minute_scale=10 means there are lines
        drawn at every 10 minute interval from start to finish
    space_between_minutes : integer (optional); default=50
        scale factor in pixel spacing between minute line markers
    colors : list (optional)
        a list of colors to choose from when coloring the nodes in the
        gantt chart

    Returns
    -------
    None
        the function does not return any value but writes out an html
        file in the same directory as the callback log path passed in

    Usage
    -----
    # import logging
    # import logging.handlers
    # from nipype.utils.profiler import log_nodes_cb

    # log_filename = 'callback.log'
    # logger = logging.getLogger('callback')
    # logger.setLevel(logging.DEBUG)
    # handler = logging.FileHandler(log_filename)
    # logger.addHandler(handler)

    # #create workflow
    # workflow = ...

    # workflow.run(plugin='MultiProc',
    #     plugin_args={'n_procs':8, 'memory':12, 'status_callback': log_nodes_cb})

    # generate_gantt_chart('callback.log', 8)
    """  # noqa: E501

    # add the html header
    html_string = """<!DOCTYPE html>
    <head>
        <style>
            #content{
                width:99%;
                height:100%;
                position:absolute;
            }

            .node{
                background-color:#7070FF;
                border-radius: 5px;
                position:absolute;
                width:20px;
                white-space:pre-wrap;
            }

            .line{
                position: absolute;
                color: #C2C2C2;
                opacity: 0.5;
                margin: 0px;
            }

            .time{
                position: absolute;
                font-size: 16px;
                color: #666666;
                margin: 0px;
            }

            .bar{
                position: absolute;
                height: 1px;
                opacity: 0.7;
            }

            .dot{
                position: absolute;
                width: 1px;
                height: 1px;
                background-color: red;
            }
            .label {
                width:20px;
                height:20px;
                opacity: 0.7;
                display: inline-block;
            }
        </style>
    </head>

    <body>
        <div id="content">
            <div style="display:inline-block;">
    """  # noqa: E501

    close_header = """
    </div>
    <div style="display:inline-block;margin-left:60px;vertical-align: top;">
        <p><span><div class="label" style="background-color:#90BBD7;"></div> Estimated Resource</span></p>
        <p><span><div class="label" style="background-color:#03969D;"></div> Actual Resource</span></p>
        <p><span><div class="label" style="background-color:#f00;"></div> Failed Node</span></p>
    </div>
    """  # noqa: E501

    # Read in json-log to get list of node dicts
    nodes_list = log_to_dict(logfile)

    # Only include nodes with timing information, and covert timestamps
    # from strings to datetimes
    try:
        nodes_list = _timing(nodes_list)
    except ProcessLookupError:
        warn(str(ProcessLookupError))
        nodes_list = []

    if not nodes_list:
        return

    for node in nodes_list:
        if "duration" not in node:
            node["duration"] = (node["finish"] - node["start"]).total_seconds()

    # Create the header of the report with useful information
    start_node = nodes_list[0]
    last_node = nodes_list[-1]
    duration = (last_node["finish"] - start_node["start"]).total_seconds()

    # Get events based dictionary of node run stats
    events = create_event_dict(start_node["start"], nodes_list)

    # Summary strings of workflow at top
    html_string += (
        "<p>Start: " + start_node["start"].strftime("%Y-%m-%d %H:%M:%S")
        + "</p>"
    )
    html_string += (
        "<p>Finish: " + last_node["finish"].strftime("%Y-%m-%d %H:%M:%S")
        + "</p>"
    )
    html_string += "<p>Duration: " + "{0:.2f}".format(duration / 60) \
                   + " minutes</p>"
    html_string += "<p>Nodes: " + str(len(nodes_list)) + "</p>"
    html_string += "<p>Cores: " + str(cores) + "</p>"
    html_string += close_header
    # Draw nipype nodes Gantt chart and runtimes
    html_string += draw_lines(
        start_node["start"], duration, minute_scale, space_between_minutes
    )
    html_string += draw_nodes(
        start_node["start"],
        nodes_list,
        cores,
        minute_scale,
        space_between_minutes,
        colors,
    )

    # Get memory timeseries
    estimated_mem_ts = calculate_resource_timeseries(
        events, "estimated_memory_gb")
    runtime_mem_ts = calculate_resource_timeseries(events, "runtime_memory_gb")
    # Plot gantt chart
    resource_offset = 120 + 30 * cores
    html_string += draw_resource_bar(
        start_node["start"],
        last_node["finish"],
        estimated_mem_ts,
        space_between_minutes,
        minute_scale,
        "#90BBD7",
        resource_offset * 2 + 120,
        "Memory",
    )
    html_string += draw_resource_bar(
        start_node["start"],
        last_node["finish"],
        runtime_mem_ts,
        space_between_minutes,
        minute_scale,
        "#03969D",
        resource_offset * 2 + 120,
        "Memory",
    )

    # Get threads timeseries
    estimated_threads_ts = calculate_resource_timeseries(
        events, "estimated_threads")
    runtime_threads_ts = calculate_resource_timeseries(
        events, "runtime_threads")
    # Plot gantt chart
    html_string += draw_resource_bar(
        start_node["start"],
        last_node["finish"],
        estimated_threads_ts,
        space_between_minutes,
        minute_scale,
        "#90BBD7",
        resource_offset,
        "Threads",
    )
    html_string += draw_resource_bar(
        start_node["start"],
        last_node["finish"],
        runtime_threads_ts,
        space_between_minutes,
        minute_scale,
        "#03969D",
        resource_offset,
        "Threads",
    )

    # finish html
    html_string += """
        </div>
    </body>"""

    # save file
    with open(logfile + ".html", "w") as html_file:
        html_file.write(html_string)


def resource_overusage_report(cblog):
    '''Function to parse through a callback log for memory and/or
    thread usage above estimates / limits.

    Parameters
    ----------
    cblog: str
        path to callback.log

    Returns
    -------
    text_report: str

    excessive: dict
    '''
    cb_dict_list = log_to_dict(cblog)
    excessive = {node['id']: [
        node['runtime_memory_gb']if node.get('runtime_memory_gb', 0)
        > node.get('estimated_memory_gb', 1) else None,
        node['estimated_memory_gb'] if node.get('runtime_memory_gb', 0)
        > node.get('estimated_memory_gb', 1) else None,
        node['runtime_threads'] - 1 if node.get('runtime_threads', 0) - 1
        > node.get('num_threads', 1) else None,
        node['num_threads'] if node.get('runtime_threads', 0) - 1
        > node.get('num_threads', 1) else None
    ] for node in [node for node in cb_dict_list if (
        node.get('runtime_memory_gb', 0) > node.get('estimated_memory_gb', 1)
        # or node.get('runtime_threads', 0) - 1 > node.get('num_threads', 1)
    )]}
    text_report = ''
    if excessive:
        text_report += 'The following nodes used excessive resources:\n'
        dotted_line = '-' * (len(text_report) - 1) + '\n'
        text_report += dotted_line
        for node in excessive:
            node_id = '\n  .'.join(node.split('.'))
            text_report += f'\n{node_id}\n'
            if excessive[node][0]:
                text_report += '      **memory_gb**\n' \
                               '        runtime > estimated\n' \
                               f'        {excessive[node][0]} ' \
                               f'> {excessive[node][1]}\n'
            # JC: I'm not convinced 'runtime_threads' and 'threads' are
            # comparable in nipype ~1.5.1
            # if excessive[node][2]:
            #     text_report += '      **threads**\n        runtime > limit\n' \
            #                    f'        {excessive[node][2]} ' \
            #                    f'> {excessive[node][3]}\n'
        text_report += dotted_line
    return text_report, excessive


def resource_report(callback_log, num_cores, logger=None):
    '''Function to attempt to warn any excessive resource usage and
    generate an interactive HTML chart.

    Parameters
    ----------
    callback_log: str
        path to callback.log

    num_cores: int
    logger: Logger
        https://docs.python.org/3/library/logging.html#logger-objects

    Returns
    -------
    None
    '''
    e_msg = ''
    try:
        txt_report = resource_overusage_report(callback_log)[0]
        if txt_report:
            with open(callback_log + '.resource_overusage.txt',
                      'w') as resource_overusage_file:
                resource_overusage_file.write(txt_report)
    except Exception as exception:  # pylint: disable=broad-except
        e_msg += f'Excessive usage report failed for {callback_log} ' \
                 f'({str(exception)})\n'
    generate_gantt_chart(callback_log, num_cores)
    if e_msg:
        if logger is not None:
            logger.warning(e_msg, exc_info=1)
        else:
            warn(e_msg, category=ResourceWarning)


def _timing(nodes_list):
    """Covert timestamps from strings to datetimes

    Parameters
    ----------
    nodes_list : list of dicts

    Returns
    -------
    nodes_list : list of dicts

    Examples
    --------
    >>> node_list = [{"start": "2022-04-16T14:43:36.286896",
    ...               "finish": "2022-04-16T14:43:36"},
    ...              {"start": "2022-04-16T14:43:38",
    ...               "finish": "2022-04-16T14:43:38.529389"}]
    >>> all(isinstance(v, str) for node in node_list for v in node.values())
    True
    >>> any(isinstance(v, datetime) for node in node_list for
    ...     v in node.values())
    False
    >>> all(isinstance(v, datetime) for node in _timing(node_list) for
    ...     v in node.values())
    True
    """
    try:
        return [_timing_timestamp(node) for node in nodes_list if
                "start" in node and "finish" in node]
    except ValueError:
        # Drop any problematic nodes
        new_node_list = []
        for node in nodes_list:
            if "start" in node and "finish" in node:
                try:
                    new_node_list.append(_timing_timestamp(node))
                except ValueError:
                    pass
        return new_node_list


def _timing_timestamp(node):
    """Convert timestamps in a node from string to datetime

    Parameters
    ----------
    node : dict

    Returns
    -------
    dict
    """
    if node is None or node.items() is None:
        raise ProcessLookupError('No logged nodes have timing information.')
    return {k: (datetime.strptime(v, "%Y-%m-%dT%H:%M:%S.%f") if '.' in v else
                datetime.fromisoformat(v)) if (k in {"start", "finish"} and
                                               isinstance(v, str)) else
            v for k, v in node.items()}

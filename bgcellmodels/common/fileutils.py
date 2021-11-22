# -*- coding: utf-8 -*-

"""
File handling utilities.
"""

import re
import json
import collections
from _ctypes import PyObj_FromPtr


def parse_json_file(filename, nonstrict=True, ordered=False):
    """
    Parse JSON file to dict.

    @param      nonstrict: bool

                Allow comments and trailing commas in json string.
    """
    with open(filename, 'r') as json_file:
        json_string = json_file.read()
        return parse_json_string(
                    json_string, nonstrict=nonstrict, ordered=ordered)


def parse_json_string(string, nonstrict=True, ordered=False):
    """
    Parse JSON string to dict.

    @param      nonstrict: bool

                Allow comments and trailing commas in json string.
    """
    if nonstrict:
        # string = validate_minify_json(string)
        string = remove_trailing_commas(remove_comments(string))
    if ordered:
        object_pairs_hook = collections.OrderedDict
    else:
        object_pairs_hook = None
    return json.loads(string, object_pairs_hook=object_pairs_hook)


def remove_comments(json_like):
    """
    Removes C-style comments from *json_like* and returns the result.

    @author     Dan McDougall <daniel.mcdougall@liftoffsoftware.com>

    Example
    -------

    >>> test_json = '''\
    {
        "foo": "bar", // This is a single-line comment
        "baz": "blah" /* Multi-line
        Comment */
    }'''
    >>> remove_comments('{"foo":"bar","baz":"blah",}')
    '{\n    "foo":"bar",\n    "baz":"blah"\n}'

    """
    comments_re = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    def replacer(match):
        s = match.group(0)
        if s[0] == '/': return ""
        return s
    return comments_re.sub(replacer, json_like)


def remove_trailing_commas(json_like):
    """
    Removes trailing commas from *json_like* and returns the result.
    
    @author     Dan McDougall <daniel.mcdougall@liftoffsoftware.com>


    Example
    -------

    >>> remove_trailing_commas('{"foo":"bar","baz":["blah",],}')
        '{"foo":"bar","baz":["blah"]}'

    """

    trailing_object_commas_re = re.compile(
        r'(,)\s*}(?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    trailing_array_commas_re = re.compile(
        r'(,)\s*\](?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    # Fix objects {} first
    objects_fixed = trailing_object_commas_re.sub("}", json_like)
    # Now fix arrays/lists [] and return the result
    return trailing_array_commas_re.sub("]", objects_fixed)


class NoIndent(object):
    """
    Wrapper for lists/dicts that should not be indented in JSON representation.
    """
    def __init__(self, value):
        self.value = value


class VariableIndentEncoder(json.JSONEncoder):
    """
    JSON encoder that lets you combine indentation with 'flat' (not indented)
    representations of objects.

    This makes the JSON file more readable for example when writing matrices
    as nested attributes.

    Credit to StackOverflow user martineau, https://stackoverflow.com/a/13252112

    @note   must write to string, not directly to file. I.e. use json.dumps
            instead of json.dump

    EXAMPLE
    -------

    >>> properties = {
    >>>     'flat_dict': NoIndent([{"x":1,"y":7}, {"x":0,"y":4}, {"x":5,"y":3}]),
    >>>     'flat_list': NoIndent([[1,2,3],[2,3,1],[3,2,1]])
    >>> }
    >>>
    >>> json.dumps(properties, cls=VariableIndentEncoder, indent=2, sort_keys=False)

    """
    FORMAT_SPEC = '@@{}@@'
    regex = re.compile(FORMAT_SPEC.format(r'(\d+)'))

    def __init__(self, **kwargs):
        # Save copy of any keyword argument values needed for use here.
        self.__sort_keys = kwargs.get('sort_keys', None)
        super(VariableIndentEncoder, self).__init__(**kwargs)


    def default(self, obj):
        # The default JSON representation of a NoIndent object is its id (memory address)
        if isinstance(obj, NoIndent):
            return self.FORMAT_SPEC.format(id(obj))
        else:
            return super(VariableIndentEncoder, self).default(obj)


    def encode(self, obj):
        format_spec = self.FORMAT_SPEC  # Local var to expedite access.
        json_repr = super(VariableIndentEncoder, self).encode(obj)  # Default JSON.

        # Replace any marked-up object ids in the JSON repr with the
        # value returned from the json.dumps() of the corresponding
        # wrapped Python object.
        for match in self.regex.finditer(json_repr):
            # Get id (memory address) of wrapper object and retrieve it
            id = int(match.group(1))
            no_indent = PyObj_FromPtr(id)
            json_obj_repr = json.dumps(no_indent.value, sort_keys=self.__sort_keys)

            # Replace the matched id string with json formatted representation
            # of the corresponding Python object.
            json_repr = json_repr.replace(
                            '"{}"'.format(format_spec.format(id)), json_obj_repr)

        return json_repr

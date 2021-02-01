#!/usr/bin/env python

import simplejson as json
import os
from datetime import datetime

class JSON(object):
    
    def __init__(self):
        """
        Construct a new JSON object.
        Creates the data dictionary.
        """
        self.data = {}
        self.filename = None

    def name(self, sampleName):
        """
        :param sampleName: The name of JSON object
        """
        date = datetime.now().strftime('%Y%m%d%H%M%S')
        self.data[sampleName] = date
        self.data['id'] = '{}_{}'.format(sampleName, date)
        self.filename = '{}.json'.format(self.data['id'])
        
    def get_id(self):
        """
            Easy getter for getting the id value
        """
        return self.data['id']
    
    def write(self, outputDir):
        """
        Write JSON to disk
    
        :param JSON: The JSON object
        :param outputDir: Full path to dictionary to write JSON object to disk
        :return out.name: returns full path of JSON object on disk
        """
        outFile = os.path.join(outputDir, self.filename)
        out = open(outFile, 'w')
        json.dump(self.data, out)
        return os.path.abspath(out.name)
    
    def add(self, dict):
        """
        Add dictionary to JSON
        
        :param dict: dictionary that can be added
        
        """
        self.data.update(dict)
        
    def pretty_print_key(self, key):
        """
        Pretty print JSON for specific key
        
        :param key: key in the JSON for printing
        
        """    
        data = json.dumps(self.data[key], sort_keys=True, indent=4 * ' ')
        print("{}\n{}".format(key, data))
        
    def pretty_print(self):
        """
        Pretty print JSON for complete JSON
        """   
        data = json.dumps(self.data, sort_keys=True, indent=4 * ' ')
        print(data)
        
    def get_dict(self):
        """
        Pretty print JSON for complete JSON
        """
        return(self.data)
        
    def open(self, json_path=""):
        """
        Open and overwrite JSON data object
        
        :param input: location of JSON text file
        
        """
        if json_path and os.path.exists(json_path):
            with open(json_path) as f:
                self.data = json.load(f)
                self.filename = os.path.basename(json_path)
    
    def add_results(self, tool, list_of_tuples):
        """
            Easy way to inject results of a tool
            into the JSON
            
            Resulting format will be:
            self.data[tool][parameter] = value
            
            :param toolname: name of the tool
            :param list_of_tuples: [(parameter,value),(parameter,value),(..),..]
        """
        # structure in JSON
        data = dict()
        data[tool] = {}
        
        for t in list_of_tuples:
            data[tool][t[0]] = t[1]
    
        self.add(data)

    def update(self, key, list_of_tuples):
        """
            Easy way to inject results of a tool
            in a update fashion into the JSON

            Resulting format will be:
            self.data[tool][parameter] = value

            :param toolname: name of the tool
            :param list_of_tuples: [(parameter,value),(parameter,value),(..),..]
        """
        # structure in JSON
        data = dict()
        data[key] = {}

        for t in list_of_tuples:
            self.data[key].update({t[0]: t[1]})

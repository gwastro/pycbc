# convenience classes for accessing sqlite tables
# the 'get_column()' method is implemented for instances of TableData
# in parallel to the existing pylal.SnglInspiralUtils functions

import sqlite3
import numpy as np

class SqliteData(object):
    """parent class for doing stuff with sqlite dbs"""
    
    def connect_db(self, db):
        self.conn = sqlite3.connect(db)
        self.curs = self.conn.cursor()

    def commit_to_db(self, db):
        self.conn.commit()

    def disconnect(self):
        self.conn.close()

class TableData(SqliteData):
    """class for accessing a single table in a db"""
    
    def __init__(self, tablename, columnlist, verbose=False):
        if len(columnlist):
            self.columns = columnlist
            self.indices = dict(zip(columnlist,range(len(columnlist))))
        else: raise ValueError("I need a column list to exist!")
        if str(tablename) == tablename:
            self.tablen = tablename
        else: raise ValueError("I need a tablename that is a string!")
        self.data = None
        self.basecommand = "SELECT "+", ".join(col for col in self.columns)+" FROM "+self.tablen
        self.command = self.basecommand
        self.restriction = ""
        self.verbose = verbose

    def restrict_rows(self, conditions=None, limit=None):
        """
        restrict the rows that will be returned from db
        this voids any previously stored data and resets 
        any previous restriction
        """
        self.data = None
        self.restriction = ""
        if len(conditions):
            whereFlag = True
            for cond in conditions:
                if str(cond) != cond: raise ValueError("conditions must be a list of strings!")
                if whereFlag:
                    self.restriction += " WHERE " # first condition needs a WHERE
                    whereFlag = False
                else:
                    self.restriction += " AND " # all subsequent conditions need an AND
                self.restriction += cond
        if limit: self.restriction += " LIMIT "+str(limit)
        if self.verbose: print self.restriction

    def retrieve_rows(self):
        """
        assigns a list of rows of db content to self.data
        """
        if self.verbose: print self.command+self.restriction
        # the result of executing is an *iterator* over rows
        self.data = list(self.curs.execute(self.command+self.restriction))
        if self.verbose: print "Shape of data is", np.shape(self.data)

    def get_column(self, colname):
        """
        returns an array of the values in a column
        either via slicing a pre-existing list or querying directly
        """
        # data may already be cached as a list
        if self.data is not None:
            colindex = self.indices[colname]
            if self.verbose: print "    Getting", colname, "data from pre-existing list"
            # allow the list to be of length 0, then return an empty array
            col = np.array([row[colindex] for row in self.data] if len(self.data)>0 else [])
        else:
            column_command = "SELECT "+colname+" FROM "+self.tablen
            if self.verbose: print column_command+self.restriction
            col = np.array(val[0] for val in self.curs.execute(column_command+self.restriction))
            if self.verbose: print "Shape of column", colname, "is", np.shape(col)
        return col


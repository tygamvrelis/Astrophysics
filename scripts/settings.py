# Global data and settings
# Author: Tyler
# Date: July 10, 2019

class Settings(object):
    """
    Globally-accessible object that aggregates parameters for this run
    """
    def __init__(self):
        self.__P_first = True
        self.__V_first = True
        self.__L_first = True
        self.__user_params_first = True
    
    @property
    def UNIT_DENSITY(self):
        return self.__P
    
    @UNIT_DENSITY.setter
    def UNIT_DENSITY(self, value):
        if not self.__P_first and value != self.__P:
            print("WARNING: Overwriting UNIT_DENSITY!")
            print("\tNew: {0}, prev: {1}".format(value, self.__P))
        else:
            self.__P_first = False
        self.__P = value
        
    @property
    def UNIT_VELOCITY(self):
        return self.__V
    
    @UNIT_VELOCITY.setter
    def UNIT_VELOCITY(self, value):
        if not self.__V_first and value != self.__V:
            print("WARNING: Overwriting UNIT_VELOCITY!")
            print("\tNew: {0}, prev: {1}".format(value, self.__V))
        else:
            self.__V_first = False
        self.__V = value

    @property
    def UNIT_LENGTH(self):
        return self.__L
    
    @UNIT_LENGTH.setter
    def UNIT_LENGTH(self, value):
        if not self.__L_first and value != self.__L:
            print("WARNING: Overwriting UNIT_LENGTH!")
            print("\tNew: {0}, prev: {1}".format(value, self.__L))
        else:
            self.__L_first = False
        self.__L = value

    @property
    def user_params(self):
        return self.__user_params
    
    @user_params.setter
    def user_params(self, params):
        if not self.__user_params_first and params != self.__user_params:
            print("WARNING: Overwriting user_params!")
            print("\tNew: {0}, prev: {1}".format(params, self.__user_params))
        else:
            self.__user_params_first = False
        self.__user_params = params

settings = Settings()
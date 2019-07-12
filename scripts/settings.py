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
        self.__T_first = True
        self.__user_params_first = True
        self.__B0_first = True
        self.__eta_first = True
    
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
    def UNIT_TIME(self):
        return self.__T
    
    @UNIT_TIME.setter
    def UNIT_TIME(self, value):
        if not self.__T_first and value != self.__T:
            print("WARNING: Overwriting UNIT_TIME!")
            print("\tNew: {0}, prev: {1}".format(value, self.__T))
        else:
            self.__T_first = False
        self.__T = value

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

    @property
    def B0(self):
        return self.__B0
    
    @B0.setter
    def B0(self, field):
        if not self.__B0_first and field != self.__B0:
            print("WARNING: Overwriting B0!")
        else:
            self.__B0_first = False
        self.__B0 = field

    @property
    def eta(self):
        return self.__eta
    
    @eta.setter
    def eta(self, field):
        if not self.__eta_first and field != self.__eta:
            print("WARNING: Overwriting eta!")
        else:
            self.__eta_first = False
        self.__eta = field

settings = Settings()
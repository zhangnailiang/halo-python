import os

def metakernelcheck():
    """
    This function creates a 'metakernel.tm' file and writes the necessary kernel
    loading commands into it.
    """
    strfirst = 'KPL/MK \n \\begindata'
    strend = (
        '\n\t\tKERNELS_TO_LOAD   = ('
        '\n\t\t\t\'./ker/de430.bsp\','
        '\n\t\t\t\'./ker/pck00010.tpc\','
        '\n\t\t\t\'./ker/naif0012.tls\','
        '\n\t\t\t\'./ker/moon_pa_de421_1900-2050.bpc\','
        '\n\t\t\t\'./ker/moon_080317.tf\','
        '\n\t\t\t\'./ker/moon_assoc_me.tf\','
        '\n\t\t\t\'./ker/earth_1962_240827_2124_combined.bpc\')'
        '\n\\begintext'
    )
    strdoc = strfirst + strend
    with open('metakernel.tm', 'w') as file:
        file.write(strdoc)
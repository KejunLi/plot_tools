#!/usr/bin/env python3
import matplotlib.colors as mc
import colorsys

# this is a color list for plot
def var_color(color, brightness_offset):
    """
    lightens the given color by multiplying (1-luminosity)
    by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    Examples:
    var_color("g", 0.3)
    var_color("#F034A3", 0.6)
    var_color((0.3,0.55,0.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    # print(colorsys.hls_to_rgb(c[0], 1-brightness_offset*(1-c[1]), c[2]))
    return(colorsys.hls_to_rgb(c[0], 1-brightness_offset*(1-c[1]), c[2]))


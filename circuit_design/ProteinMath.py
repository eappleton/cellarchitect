import matplotlib
import matplotlib.pyplot as plt
import dnaplotlib as dpl
from matplotlib import gridspec
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch
from matplotlib.path import Path
import matplotlib.path as mpath
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches
import math
import random
import csv
import pdb

base_colors = {"\N{GREEK CAPITAL LETTER PHI}"+'C31': (0.694, 0.49, .772),
               'Bxb1': (.152, 0.69, 0.643),
               'TP901': (0.192, 0.866, 0.192),
               'iRFP670': (0.768, 0.109, .109),
               'EBFP2': (.109, 0.47, 0.768),
               'mEmerald': (.184, 0.505, 0.113),
               'mScarlet-I': (.937, .101, .101),
               "\N{GREEK CAPITAL LETTER PHI}"+'BT1': (.937, .101, 0.384),
               'int12': (.929, 0.47, 0.25),
               "\N{GREEK CAPITAL LETTER PHI}"+'Joe': (0.25, 0.929, 0.607),
               'sprA': (0.921, 0.913, 0),
               "\N{GREEK CAPITAL LETTER PHI}"+'Rv1':(1,.2,.2),
               'int11':(.2,1,.2),
               'int8':(.2,.2,1),
               'int4':(.6,.9,0),
               'int9':(.1,.3,.5),
               'l': (0.952, 0.584, 0.086),
               'm': (0.984, 0.035, 0.917),
               'gene': (0.850, 0.788, 0.788),
               'unknown': (0.850, 0.788, 0.788)}
#background=(.9137,.9137,.9137)
background=(1,1,1)

letter_to_color = {'a': 'khaki',
               'b': 'red',
               'c': 'light green',
               'd': 'light blue',
               'e': 'yellow',
               'f': 'violet',
               'g': 'dark blue',
               'h': 'cyan',
               'i': 'dark green',
               'j': 'black',
               'k': 'brown',
               'l': 'orange',
               'm': 'pink',
               'gene': 'grey',
               'numb': 'gold',
               'unknown': 'white'}

recSiteLookupBP={"1": 'RecombinaseSite',
               "3": 'RecombinaseSite_round',
               "4": 'RecombinaseSite_arrow',
               "2": 'RecombinaseSite_flying',
               "5": 'RecombinaseSite_block',
               "6": 'RecombinaseSite_nepalese_flag'}

recSiteLookupLR={"1": 'RecombinaseSite2',
               "3": 'RecombinaseSite_round2',
               "4": 'RecombinaseSite_arrow2',
               "2": 'RecombinaseSite_flying2',
               "5": 'RecombinaseSite_block2',
               "6": 'RecombinaseSite_nepalese_flag2'}
def makeCircuits():
    x=["ABC","BC","AC","AB","A","B","C",""]
    count=0
    for i in range(8):
        for j in range(i,8):
            f=open("Protein_Math/Circuit"+str(count)+".csv","w")
            f.write("pEf1alpha_promoter_F\nPhiC31_1_attB_site_F\n")
            for s in x[i]:
                f.write(s+"_cassette_F\n")
            f.write("sv40_terminator_F\nPhiC31_1_attP_site_F\n")
            for s in x[j]:
                f.write(s+"_cassette_F\n")
            f.write("sv40_terminator_F")
            f.close()
            count+=1


def sbol_recombinase1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start, y_lower), 
                  (start, y_upper),
                  (end,0)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end
    
    
    
def sbol_recombinase2 (ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - LS mode
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
        else:
            if 'color' in list(opts.keys()):
                r2 = float(color[0]) / 2
                g2 = float(color[1]) / 2
                b2 = float(color[2]) / 2
                color2 = (r2,g2,b2)
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp = color
        color = color2
        color2 = temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start, y_lower), 
                 (start, y_upper),
                  (end,0)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")]) 
    midpoint = (end + start) / 2
    hypotenuse = math.sqrt( (y_extent/2)**2 + (x_extent)**2 )
    hypotenuse2 = hypotenuse / 2
    cosineA = (y_extent/2) / hypotenuse
    f = hypotenuse2 * cosineA
    p2 = Polygon([(midpoint, -1*f), 
                  (midpoint, f),
                  (end,0)],
                  edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=12, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    ax.add_patch(p2)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_round1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(start, y_lower)))
    p1.append((Path.LINETO,(start, y_upper)))
    #p1.append((Path.CURVE3,(2*end-start,0)))
    #p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.CURVE4, ((end-.25*start)/.75,y_upper)))
    p1.append((Path.CURVE4, ((end-.25*start)/.75,y_lower)))
    p1.append((Path.CURVE4, (start,y_lower)))
    p1.append((Path.LINETO,(start, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end
    
def sbol_recombinase_round2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp=color
        color = color2
        color2=temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(start, y_lower)))
    p1.append((Path.LINETO,(start, y_upper)))
    #p1.append((Path.CURVE3,(2*end-start,0)))
    #p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.CURVE4, ((end-.25*start)/.75,y_upper)))
    p1.append((Path.CURVE4, ((end-.25*start)/.75,y_lower)))
    p1.append((Path.CURVE4, (start,y_lower)))
    p1.append((Path.LINETO,(start, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    midT=.1464466094
    y_upper2=(1-midT)**3*y_upper+3*(1-midT)**2*midT*y_upper+3*(1-midT)*midT**2*y_lower+midT**3*y_lower
    y_upper3=.6*y_upper2
    y_lower2=-y_upper2
    midpoint=(start+end)/2
    p1=[]
    p1.append((Path.MOVETO,(midpoint, y_lower2)))
    p1.append((Path.LINETO,(midpoint, y_upper2)))
    #p1.append((Path.CURVE3,(2*end-start,0)))
    #p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.CURVE4, ((end-.25*midpoint)/.75,y_upper3)))
    p1.append((Path.CURVE4, ((end-.25*midpoint)/.75,-y_upper3)))
    p1.append((Path.CURVE4, (midpoint,y_lower2)))
    p1.append((Path.LINETO,(midpoint, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_arrow1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  (end,0), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end
    
def sbol_recombinase_arrow2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp=color
        color = color2
        color2=temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  (end,0), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  ((start+end)/2,y_upper/2),((start+end)/2,y_lower/2), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_pinch1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(start, y_lower)))
    p1.append((Path.LINETO,(start, y_upper)))
    p1.append((Path.CURVE3,(start+(end-start)/3,0)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,(start+(end-start)/3,0)))
    p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.LINETO,(start, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_pinch2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp=color
        color = color2
        color2=temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    inflection=start+(end-start)/3
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(start, y_lower)))
    p1.append((Path.LINETO,(start, y_upper)))
    p1.append((Path.CURVE3,(inflection,0)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,(inflection,0)))
    p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.LINETO,(start, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    midT=2.5**.5-1
    y_upper2=(1-midT)**2*y_upper
    midpoint=(start+end)/2
    inflection2x=midpoint+(end-midpoint)/3
    print(start)
    print(inflection2x)
    p1=[]
    p1.append((Path.MOVETO,(midpoint, -y_upper2)))
    p1.append((Path.LINETO,(midpoint, y_upper2)))
    p1.append((Path.CURVE3,(inflection2x,0)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,(inflection2x,0)))
    p1.append((Path.CURVE3,(midpoint, -y_upper2)))
    p1.append((Path.LINETO,(midpoint, 0)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_flying1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(end,0)))
    p1.append((Path.CURVE3,((start+end)/2,y_lower)))
    p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.CURVE3,((start+end)/2,0)))
    p1.append((Path.CURVE3,(start, y_upper)))
    p1.append((Path.CURVE3,((start+end)/2,y_upper)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,((start+end)/2,y_lower)))
    p1.append((Path.CURVE3,(start, y_lower)))
    #p1.append((Path.CURVE3,((start+end)/2,0)))
    #p1.append((Path.CURVE3,(start, y_upper)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_flying2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp=color
        color = color2
        color2=temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1=[]
    p1.append((Path.MOVETO,(end,0)))
    p1.append((Path.CURVE3,((start+end)/2,y_lower)))
    p1.append((Path.CURVE3,(start, y_lower)))
    p1.append((Path.CURVE3,((start+end)/2,0)))
    p1.append((Path.CURVE3,(start, y_upper)))
    p1.append((Path.CURVE3,((start+end)/2,y_upper)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,((start+end)/2,y_lower)))
    p1.append((Path.CURVE3,(start, y_lower)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    midpoint=(start+end)/2
    y_upper2=.75*y_upper
    p1=[]
    p1.append((Path.MOVETO,(end,0)))
    p1.append((Path.CURVE3,((midpoint+end)/2,-.7*y_upper2)))
    p1.append((Path.CURVE3,(midpoint, -y_upper2)))
    p1.append((Path.LINETO,(midpoint,y_upper2)))
    p1.append((Path.CURVE3,((midpoint+end)/2,.7*y_upper2)))
    p1.append((Path.CURVE3,(end, 0)))
    p1.append((Path.CURVE3,((midpoint+end)/2,-.7*y_upper2)))
    p1.append((Path.CURVE3,(midpoint, -y_upper2)))
    codes, verts= zip(*p1)
    path=mpath.Path(verts, codes)
    patch=patches.PathPatch(path, edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(patch)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_block1(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        color = color2
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  (start+(end-start)*.6, y_upper), (end,0), (start+(end-start)*.6, y_lower), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end
    
def sbol_recombinase_block2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - forward direction
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp=color
        color = color2
        color2=temp
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
    # Draw the site
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  (start+(end-start)*.6, y_upper), (end,0), (start+(end-start)*.6, y_lower), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=11, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    p1 = Polygon([(start+(end-start)*.4, 0),  (start, y_upper),
                  ((start+end)/2, y_upper), ((start+end)/2, y_lower), (start, y_lower)],
                  edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=12, 
                  path_effects=[Stroke(joinstyle="miter")])
    ax.add_patch(p1)
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end
    
def sbol_recombinase_nepalese_flag(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - alternative design, nepalese flag
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
                
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        middle_point = end + (start-end)*0.8
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
        middle_point = start + (end-start)*0.2
                
    ## Draw the site    
    verts = [(start, y_lower), # left, bottom
             (start, y_upper), # left, top
             (end, (y_upper)*0.7), # right, top
             (middle_point, (y_lower+y_upper)*0.5), # middle cut interior point
             (end, (y_lower)*0.7), # right, bottom
             (0., 0.) # ignored
            ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=linewidth, zorder=11)
    ax.add_patch(patch)

    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end

def sbol_recombinase_nepalese_flag2(ax, type, num, start, end, prev_end, scale, linewidth, opts):
    """ 
    SBOL recombinase site renderer - alternative design, nepalese flag LR mode
    """
    # Default parameters
    color = (0,0,0)
    color2 = (0,0,0)
    start_pad = 0.0
    end_pad = 0.0
    x_extent = 6.0
    y_extent = 6.0
    linestyle = '-'
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'start_pad' in list(opts.keys()):
            start_pad = opts['start_pad']
        if 'end_pad' in list(opts.keys()):
            end_pad = opts['end_pad']
        if 'x_extent' in list(opts.keys()):
            x_extent = opts['x_extent']
        if 'y_extent' in list(opts.keys()):
            y_extent = opts['y_extent']
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'color2' in list(opts.keys()):
            color2 = opts['color2']
                
    # Check direction add start padding
    final_end = end
    final_start = prev_end
    y_lower = -1 * y_extent/2
    y_upper = y_extent/2
    if start > end:
        start = prev_end+end_pad+x_extent+linewidth
        end = prev_end+end_pad
        final_end = start+start_pad
        temp = color
        color = color2
        color2 = temp
        middle_point = end + (start-end)*0.8
    else:
        start = prev_end+start_pad+linewidth
        end = start+x_extent
        final_end = end+end_pad
        middle_point = start + (end-start)*0.3
        
    ## Draw the site    
    verts = [(start, y_lower), # left, bottom
             (start, y_upper), # left, top
             (end, (y_upper)*0.7), # right, top
             (middle_point, (y_lower+y_upper)*0.5), # middle cut interior point
             (end, (y_lower)*0.7), # right, bottom
             (0., 0.) # ignored
            ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=linewidth, zorder=11)
    ax.add_patch(patch)
    
    verts = [(start, y_lower), # left, bottom
             (start, y_upper), # left, top
             ((start+end)/2, (y_upper)*0.85), # right, top
             ((start+end)/2, .35*y_upper), # middle cut interior point
             (middle_point, 0),
             ((start+end)/2, -.35*y_upper),
             ((start+end)/2, -(y_upper)*0.85),
             (0., 0.) # ignored
            ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color2, lw=linewidth, zorder=11)
    ax.add_patch(patch)
    
    # Add a label if needed
    if opts != None and 'label' in list(opts.keys()):
        if final_start > final_end:
            write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
        else:
            write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
    # Return the final start and end positions to the DNA renderer
    if final_start > final_end:
        return prev_end, final_start
    else:
        return prev_end, final_end


def write_label (ax, label_text, x_pos, opts=None):
    """ Renders labels on parts.
    """
    zorder_add = 0.0
    y_offset = 0.0
    label_style = 'normal'
    label_size = 7
    label_y_offset = 0
    label_x_offset = 0
    label_color = (0,0,0)
    label_rotation = 0
    if opts != None:
        if 'zorder_add' in list(opts.keys()):
            zorder_add = opts['zorder_add']
        if 'y_offset' in list(opts.keys()):
            y_offset = opts['y_offset']
        if 'label_style' in list(opts.keys()):
            label_style = opts['label_style']
        if 'label_size' in list(opts.keys()):
            label_size = opts['label_size']
        if 'label_y_offset' in list(opts.keys()):
            label_y_offset = opts['label_y_offset']
        if 'label_x_offset' in list(opts.keys()):
            label_x_offset = opts['label_x_offset']
        if 'label_color' in list(opts.keys()):
            label_color = opts['label_color']
        if 'label_rotation' in list(opts.keys()):
            label_rotation = opts['label_rotation']
    ax.text(x_pos+label_x_offset, label_y_offset+y_offset, label_text, horizontalalignment='center',
            verticalalignment='center', fontsize=label_size, fontstyle=label_style, 
            color=label_color, rotation=label_rotation, zorder=30+zorder_add)
    
def flip_arrow(ax, type, num, from_part, to_part, scale, linewidth, arc_height_index, opts):
    """ Regulation arcs for recombinase sites
    """
    # Default parameters
    color = (0.0,0.0,0.0)
    arcHeightStart = 10
    arcHeightEnd = 10
    linewidth = 2.0
    
    # Update default parameters if provided
    if opts != None:
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'arc_height_start' in list(opts.keys()):
            arcHeightStart = opts['arc_height_start']
        if 'arc_height_end' in list(opts.keys()):
            arcHeightEnd = opts['arc_height_end']
    start = (from_part['start'] + from_part['end']) / 2
    end   = (to_part['start']   + to_part['end']) / 2
    # Check direction and draw arc
    if start > end:
        arcHeightStart = -arcHeightStart
        arcHeightEnd = -arcHeightEnd
    ax.annotate('', (end, arcHeightEnd), (start, arcHeightStart), 
                ha="right", va="center", size=8, 
                arrowprops=dict(arrowstyle='->',connectionstyle="arc3,rad=-.4",lw=linewidth, color=color))
    
    
def dark (col, fac=2.0):
    return (col[0]/fac, col[1]/fac, col[2]/fac)

def greek(name):
    if "Phi" in name:
        return name[:name.index("Phi")]+"\N{GREEK CAPITAL LETTER PHI}"+name[name.index("Phi")+3:]
    if "alpha" in name:
        return name[:name.index("alpha")]+"\N{GREEK SMALL LETTER ALPHA}"+name[name.index("alpha")+5:]
    return None

def getColor(name):
    if name not in base_colors:
        if "Numb" in name:
            for x in base_colors:
                if x in name:
                    base_colors[name]=base_colors[x]
                    break
        if name not in base_colors:
            base_colors[name]=generateColors(100, 2.5)[0]
    return base_colors[name]

def grn_dnaplotlib_mapping(seq, highlight_part=None):
    
    seq_to_plot = []
    geneLabelYOffset=-.5
    textWidthFactor=3
    textWidthAdjustment=8
    textSize=10
    vertOffset=17
    for element in seq:
        first=element.find("_")
        second=element.find("_",first+1)
        third=element.find("_",second+1)
        name=element[:first]
        name2=greek(name)
        if name2:
            name=name2
        siteID=element[first+1:second]
        siteType=element[second+1:third]
        if 'promoter_F' in element:
            seq_to_plot.append({'type':'Promoter', 'name':'prom', 'fwd':True, "opts":{'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':vertOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4}})
        elif 'promoter_R' in element:
            seq_to_plot.append({'type':'Promoter', 'name':'prom', 'fwd':False, "opts":{'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':-vertOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4}})
        elif 'terminator_F' in element :
            seq_to_plot.append({'type':'Terminator', 'name':'term', 'fwd':True,"opts":{'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':vertOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4}})
        elif 'terminator_R' in element:
            seq_to_plot.append({'type':'Terminator', 'name':'term', 'fwd':False,"opts":{'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':-vertOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4}})
        elif 'cassette_F' in element:
            if highlight_part:
                if element in highlight_part:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                        'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
                    
                else:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                        'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
            else:        
                seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                    'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name, 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
        elif 'cassette_R' in element:
            if highlight_part:
                if element in highlight_part:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                'fwd':False, 'opts':{'color':getColor(name), 
                                                    'label':name, 
                                                    'label_x_offset':0, 'label_y_offset':geneLabelYOffset, 
                                                    'label_style':'italic',
                                                    'label_size':textSize,
                                                    'linewidth':4,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
                else:       
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                        'fwd':False, 'opts':{'color':getColor(name), 
                                                    'label':name, 
                                                    'label_x_offset':0, 'label_y_offset':geneLabelYOffset, 
                                                    'label_style':'italic',
                                                    'label_size':textSize,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
            else:
                seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                        'fwd':False, 'opts':{'color':getColor(name), 
                                                    'label':name, 
                                                    'label_x_offset':0, 'label_y_offset':geneLabelYOffset, 
                                                    'label_style':'italic',
                                                    'label_size':textSize,'x_extent':textWidthFactor*len(name)+textWidthAdjustment}})
        elif 'RDF_F' in element:
            if highlight_part:
                if element in highlight_part:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
                else:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                    'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
            else:        
                seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                    'fwd':True, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 'label_rotation':0, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
        elif 'RDF_R' in element:
            if highlight_part:
                if element in highlight_part:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                'fwd':False, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,
                                                     'linewidth':4,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
                else:
                    seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                    'fwd':False, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
            else:
                seq_to_plot.append({'type':'CDS', 'name':'cds', 
                                    'fwd':False, 'opts':{'color':getColor(name), 
                                                     'label':name+' RDF', 'label_x_offset':0, 
                                                     'label_y_offset':geneLabelYOffset, 
                                                     'label_style':'italic',
                                                     'label_size':textSize,'x_extent':textWidthFactor*(len(name)+4)+textWidthAdjustment}})
        elif 'attB_site_F' in element or 'attP_site_F' in element:
            seq_to_plot.append({'type':recSiteLookupBP[siteID],  'name':'a1',  
                                 'fwd':True,'opts':{'color':getColor(name), 
                                           'color2':getColor(name), 
                                           'x_extent':16, 'y_extent':12, 
                                           'start_pad':3, 'end_pad':3,
                                           'label':name+"\n"+siteType+" "+siteID, 'label_x_offset':0, 
                                           'label_y_offset':-vertOffset, 
                                           'label_style':'italic',
                                           'label_size':textSize
                                                   }})
        elif 'attB_site_R' in element or 'attP_site_R' in element:
            seq_to_plot.append({'type':recSiteLookupBP[siteID],  'name':'a1',  
                                'fwd':False,'opts':{'color':getColor(name), 
                                                      'color2':getColor(name), 
                                                      'x_extent':16, 'y_extent':12, 
                                                      'start_pad':3, 'end_pad':3,
                                                      'label':name+"\n"+siteType+" "+siteID, 'label_x_offset':0, 
                                                      'label_y_offset':-vertOffset, 
                                                      'label_style':'italic',
                                                      'label_size':textSize
                                                   }})
        elif 'attL_site_F' in element or 'attR_site_F' in element:
            seq_to_plot.append({'type':recSiteLookupLR[siteID],  'name':'a1f',  
                                'fwd':True,   'opts':{'color':getColor(name), 
                                                      'color2':dark(getColor(name)), 
                                                      'x_extent':16, 'y_extent':12, 
                                                      'start_pad':3, 'end_pad':3,
                                                      'label':name+"\n"+siteType+" "+siteID, 'label_x_offset':0, 
                                                      'label_y_offset':-vertOffset, 
                                                      'label_style':'italic',
                                                      'label_size':textSize
                                                     }})
        elif 'attL_site_R' in element or 'attR_site_R' in element:
            seq_to_plot.append({'type':recSiteLookupLR[siteID],  'name':'a1f',  
                                'fwd':False,  'opts':{'color':getColor(name), 
                                                      'color2':dark(getColor(name)), 
                                                      'x_extent':16, 'y_extent':12, 
                                                      'start_pad':3, 'end_pad':3,
                                                      'label':name+"\n"+siteType+" "+siteID, 'label_x_offset':0, 
                                                      'label_y_offset':-vertOffset, 
                                                      'label_style':'italic',
                                                      'label_size':textSize
                                                     }}) 
    return seq_to_plot


def opposite_site(site):
    
    if site == 'BP':
        return 'LR'
    elif site == 'LR':
        return 'BP'
    else:
        raise ValueError('site must be BP or LR in string format')
        
def drawStates(configs, edges, outFile):
    #print(configs)
    for config in configs:
        for x in config:
            count=len(x)-1
            while count>-1:
                if x[count]=="":
                    del x[count]
                count+=-1
        
    # Create the DNAplotlib renderer
    dr = dpl.DNARenderer(linewidth=2.0)
    # Use default renderers and append our custom ones for recombinases
    reg_renderers = dr.std_reg_renderers()
    reg_renderers['Connection'] = flip_arrow
    part_renderers = dr.SBOL_part_renderers()
    part_renderers['RecombinaseSite'] = sbol_recombinase1
    part_renderers['RecombinaseSite2'] = sbol_recombinase2
    part_renderers['RecombinaseSite_round']= sbol_recombinase_round1
    part_renderers['RecombinaseSite_round2']= sbol_recombinase_round2
    part_renderers['RecombinaseSite_arrow']= sbol_recombinase_arrow1
    part_renderers['RecombinaseSite_arrow2']= sbol_recombinase_arrow2
    part_renderers['RecombinaseSite_flying']= sbol_recombinase_flying1
    part_renderers['RecombinaseSite_flying2']= sbol_recombinase_flying2
    part_renderers['RecombinaseSite_block']= sbol_recombinase_block1
    part_renderers['RecombinaseSite_block2']= sbol_recombinase_block2
    part_renderers['RecombinaseSite_nepalese_flag'] = sbol_recombinase_nepalese_flag
    part_renderers['RecombinaseSite_nepalese_flag2'] = sbol_recombinase_nepalese_flag2
    
    if len(configs)<12:
        numLines=0
        spacers=1
        for config in configs:
            for x in config:
                numLines+=len(x)//15
                if len(x)%15>0:
                    numLines+=1
            spacers+=1
        
        width=15
        #print(numLines)
        height=(2.5*numLines)+(0.5*spacers)
        fig = plt.figure(figsize=(width,height), facecolor=background)
        gs = gridspec.GridSpec(5*numLines+spacers, 6)
        count=1
        coordinates=[]
        for config in configs:
            startCount=count
            for x in config:
                spot=0
                while spot<len(x):
                    design = grn_dnaplotlib_mapping(x[spot:min(spot+15,len(x))])
                    ax_dna = plt.subplot(gs[count:count+5,1:-1])
                    start, end = dr.renderDNA(ax_dna, design, part_renderers, regs=None, reg_renderers=None)
                    ax_dna.set_xlim([start, end])
                    ax_dna.set_ylim([-30,30])
                    ax_dna.set_aspect('equal')
                    ax_dna.set_xticks([])
                    ax_dna.set_yticks([])
                    ax_dna.axis('off')
                    count+=5
                    spot+=15
            coordinates.append((startCount+count)/2)
            count+=1
        left=plt.subplot(gs[:,0])
        right=plt.subplot(gs[:,5])
        left.set_ylim([5*numLines+spacers,0])
        right.set_ylim([5*numLines+spacers,0])
        left.set_xlim([(5*numLines+spacers)*(width/height)/6,0])
        right.set_xlim([0,(5*numLines+spacers)*(width/height)/6])
        left.set_xticks([])
        left.set_yticks([])
        left.axis('off')
        right.set_xticks([])
        right.set_yticks([])
        right.axis('off')
        drawEdges2(edges, left, right, coordinates)
        #fig.savefig(outFile, transparent=False)
        #fig.savefig(outFile[:-3]+"ps", papertype="letter",format="ps", transparent=False)
        fig.savefig(outFile[:-3]+'png', dpi=300, facecolor=background)
    else:
        pagenum=0
        for i in range(len(configs)):
            numLines=0
            spacers=1
            for x in configs[i]:
                numLines+=len(x)//15
                if len(x)%15>0:
                    numLines+=1
            width=15
            height=(2.5*numLines)+(0.5*spacers)
            fig = plt.figure(figsize=(width,height))
            gs = gridspec.GridSpec(5*numLines+spacers, 6)
            count=1
            for x in configs[i]:
                spot=0
                while spot<len(x):
                    design = grn_dnaplotlib_mapping(x[spot:min(spot+15,len(x))])
                    ax_dna = plt.subplot(gs[count:count+5,1:-1])
                    start, end = dr.renderDNA(ax_dna, design, part_renderers, regs=None, reg_renderers=None)
                    ax_dna.set_xlim([start, end])
                    ax_dna.set_ylim([-30,30])
                    ax_dna.set_aspect('equal')
                    ax_dna.set_xticks([])
                    ax_dna.set_yticks([])
                    ax_dna.axis('off')
                    count+=5
                    spot+=15
            left=plt.subplot(gs[:,0])
            right=plt.subplot(gs[:,5])
            left.set_ylim([5*numLines+spacers,0])
            right.set_ylim([5*numLines+spacers,0])
            left.set_xlim([(5*numLines+spacers)*(width/height)/6,0])
            right.set_xlim([0,(5*numLines+spacers)*(width/height)/6])
            left.set_xticks([])
            left.set_yticks([])
            left.axis('off')
            right.set_xticks([])
            right.set_yticks([])
            right.axis('off')
            fig.savefig(outFile[:-4]+"_"+str(pagenum)+'.png', dpi=300)
            pagenum+=1
        
def drawEdges(edges, left, right, coordinates):
    width=left.get_xlim()[0]
    colors=generateColors(len(edges))
    for i,edge in enumerate(edges):
        if edge[0]<edge[1]:
            midpoint=(((edge[1]-edge[0])/len(coordinates))*width, (coordinates[edge[1]]+coordinates[edge[0]])/2)
            pathRight= []
            pathRight.append((Path.MOVETO,(0,coordinates[edge[0]]+.5)))
            pathRight.append((Path.CURVE3,midpoint))
            pathRight.append((Path.CURVE3,(0, coordinates[edge[1]]-.5)))
            codes, verts = zip(*pathRight)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            right.add_patch(patch)
            makeArrow(0,coordinates[edge[1]]-.5, math.atan((coordinates[edge[1]]-.5-midpoint[1])/(.87*midpoint[0])), right, colors[i])
            right.text(midpoint[0], midpoint[1], cleanProb(edge[2]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.7, ec="None"))
        elif edge[0]>edge[1]:
            midpoint=(((edge[0]-edge[1])/len(coordinates))*width, (coordinates[edge[1]]+coordinates[edge[0]])/2)
            pathLeft= []
            pathLeft.append((Path.MOVETO,(0,coordinates[edge[0]]-.5)))
            pathLeft.append((Path.CURVE3,midpoint))
            pathLeft.append((Path.CURVE3,(0, coordinates[edge[1]]+.5)))
            codes, verts = zip(*pathLeft)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            left.add_patch(patch)
            makeArrow(0,coordinates[edge[1]]+.5, math.atan((coordinates[edge[1]]+.5-midpoint[1])/(.87*midpoint[0])), left, colors[i])
            left.text(midpoint[0], midpoint[1], cleanProb(edge[2]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.7, ec="None"))
        else:
            pathLeft= []
            midpoint=(width/2,coordinates[edge[0]]-2)
            midpoint0=(width/2,coordinates[edge[0]]+2)
            pathLeft.append((Path.MOVETO,(0,coordinates[edge[0]]+.2)))
            pathLeft.append((Path.CURVE4,midpoint0))
            pathLeft.append((Path.CURVE4,midpoint))
            pathLeft.append((Path.CURVE4,(0,coordinates[edge[0]]-.2)))
            codes, verts = zip(*pathLeft)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            left.add_patch(patch)
            makeArrow(0,coordinates[edge[1]]-.2, math.atan((coordinates[edge[1]]-.2-midpoint[1])/(1.05*midpoint[0])), left, colors[i])
            left.text(midpoint[0], (midpoint[1]+midpoint0[1])/2+.5, cleanProb(edge[2]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.7, ec="None"))

def drawEdges2(edges,left, right, coordinates):
    ySpan=max(coordinates)-min(coordinates)
    arrows=[]
    #pdb.set_trace()
    width=left.get_xlim()[0]
    colors=generateColors(len(edges))
    breakys=[]
    breakxs=[0.0]*len(edges)
    points=[]
    for i,edge in enumerate(edges):
        mid1=(coordinates[edge[1]]+coordinates[edge[0]])/2
        mid2=(coordinates[edge[2]]+coordinates[edge[0]])/2
        if abs(mid1-coordinates[edge[0]])>abs(mid2-coordinates[edge[0]]):
            breaky=mid1
        elif abs(mid1-coordinates[edge[0]])==abs(mid2-coordinates[edge[0]]):
            breaky=max(mid1,mid2)
        else:
            breaky=mid2
        breakys.append(breaky)
        breakx=(abs(mid1-coordinates[edge[0]])+abs(mid2-coordinates[edge[0]]))/(ySpan)*width*.35+.35*width
        if (edge[1]<edge[0] and edge[2]>edge[0]) or (edge[1]>edge[0] and edge[2]<edge[0]):
            adjustment=1
        else:
            adjustment=0
        points.append([breaky,adjustment,breakx,i])
    points.sort()
    spot=0
    while spot<len(points):
        increment=1
        while spot+increment<len(points) and points[spot][0]==points[spot+increment][0]:
            increment+=1
        if increment>1:
            offset=.35*width/(increment-1)
            newx=.35*width
            for i in range(spot, spot+increment):
                points[i][2]=newx
                newx+=offset
        spot+=increment
    for point in points:
        breakxs[point[3]]=point[2]
    for i, edge in enumerate(edges):
        #pdb.set_trace()
        breaky=breakys[i]
        breakx=breakxs[i]
        if edge[0]==edge[1]==edge[2]:
            p0=(0, coordinates[edge[0]])
            p1=(0.6*width, p0[1]-0.4*width)
            p2=(0.6*width, p0[1])
            left.text(0.2*width, p0[1]-0.15*width, cleanProb(edge[3]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor=background, alpha=0.7, ec="None"))
            pathLeft=[]
            pathLeft.append((Path.MOVETO, p0))
            pathLeft.append((Path.CURVE3,p1))
            pathLeft.append((Path.CURVE3,p2))
            a1=(0.6*width, p0[1]+0.4*width)
            a2=(0, p0[1]+0.5)
            t=approximate(p2, a1, a2, 0.2*width)
            dx,dy=dydx(p2, a1, a2, t)
            a3=quad(p2,a1,a2,t)
            deltaX=0.6*width-a3[0]
            deltaY=deltaX*(dy/dx)
            a4=(0.6*width, a3[1]+deltaY)
            pathLeft.append((Path.CURVE3,a4))
            pathLeft.append((Path.CURVE3,a3))
            arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
            a1=(0.6*width, p0[1]+0.6*width)
            a2=(0, p0[1]+0.5)
            pathLeft.append((Path.MOVETO,p2))
            t=approximate(p2, a1, a2, 0.2*width)
            dx,dy=dydx(p2, a1, a2, t)
            a3=quad(p2,a1,a2,t)
            deltaX=0.6*width-a3[0]
            deltaY=deltaX*(dy/dx)
            a4=(0.6*width, a3[1]+deltaY)
            pathLeft.append((Path.CURVE3,a4))
            pathLeft.append((Path.CURVE3,a3))
            arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
            codes, verts = zip(*pathLeft)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            left.add_patch(patch)
        elif breaky>coordinates[edge[0]]:
            p0=(0, coordinates[edge[0]])
            p1=(breakx, (breaky+coordinates[edge[0]])/2)
            p2=(breakx,breaky)
            textPoint=quad(p0,p1,p2,breakx/width)
            right.text(textPoint[0]+.1*width, textPoint[1], cleanProb(edge[3]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor=background, alpha=0.7, ec="None"))
            pathRight=[]
            pathRight.append((Path.MOVETO, p0))
            pathRight.append((Path.CURVE3,p1))
            pathRight.append((Path.CURVE3,p2))
            if edge[1]==edge[2]:
                a1=(breakx,breaky+(coordinates[edge[1]]-breaky)/2)
                a2=(0, coordinates[edge[1]]-.5)
                t=approximate(p2, a1, a2, 0.2*width)
                dx,dy=dydx(p2, a1, a2, t)
                a3=quad(p2,a1,a2,t)
                deltaX=breakx-a3[0]
                deltaY=deltaX*(dy/dx)
                a4=(breakx, a3[1]+deltaY)
                pathRight.append((Path.CURVE3,a4))
                pathRight.append((Path.CURVE3,a3))
                arrows.append((a3[0], a3[1], math.atan(-dy/dx), right, colors[i]))
                a1=(breakx,coordinates[edge[1]]-.5)
                a2=(0, coordinates[edge[1]]-.5)
                t=approximate(p2, a1, a2, 0.2*width)
                dx,dy=dydx(p2, a1, a2, t)
                a3=quad(p2,a1,a2,t)
                deltaX=breakx-a3[0]
                deltaY=deltaX*(dy/dx)
                a4=(breakx, a3[1]+deltaY)
                pathRight.append((Path.MOVETO,p2))
                pathRight.append((Path.CURVE3,a4))
                pathRight.append((Path.CURVE3,a3))
                arrows.append((a3[0], a3[1], math.atan(-dy/dx), right, colors[i]))
            else:
                if coordinates[edge[1]]>breaky:
                    a1=(breakx,(breaky+coordinates[edge[1]])/2)
                    a2=(0, coordinates[edge[1]]-.5)
                    t=approximate(p2, a1, a2, 0.2*width)
                    dx,dy=dydx(p2, a1, a2, t)
                    a3=quad(p2,a1,a2,t)
                    deltaX=breakx-a3[0]
                    deltaY=deltaX*(dy/dx)
                    a4=(breakx, a3[1]+deltaY)
                    pathRight.append((Path.CURVE3,a4))
                    pathRight.append((Path.CURVE3,a3))
                    arrows.append((a3[0], a3[1], math.atan(-dy/dx), right, colors[i]))
                else:
                    if edge[1]>=edge[0]:
                        if coordinates[edge[1]]==breaky:
                            a1=(breakx, breaky+breakx*.75)
                            a2=(breakx/2, breaky+breakx*.75)
                            a3=(breakx/2, breaky+breakx*.25) 
                        else:
                            a1=(breakx, breaky+breakx/2)
                            a2=(breakx/2, breaky+breakx/2)
                            a3=(breakx/2, breaky)
                        a4=(breakx/2, (breaky+coordinates[edge[1]])/2)
                    else:
                        a1=(breakx, breaky+(width-breakx)/2)
                        a2=((width+breakx)/2, breaky+(width-breakx)/2)
                        a3=((width+breakx)/2, breaky)
                        a4=((width+breakx)/2, (breaky+coordinates[edge[1]])/2)
                    a5=(0, coordinates[edge[1]]+.5)
                    pathRight.append((Path.CURVE4, a1))
                    pathRight.append((Path.CURVE4, a2))
                    pathRight.append((Path.CURVE4, a3))
                    t=approximate(a3, a4, a5, 0.2*width)
                    dx,dy=dydx(a3, a4, a5, t)
                    a6=quad(a3,a4,a5,t)
                    deltaX=a3[0]-a6[0]
                    deltaY=deltaX*(dy/dx)
                    a7=(a3[0], a6[1]+deltaY)
                    pathRight.append((Path.CURVE3,a7))
                    pathRight.append((Path.CURVE3,a6))
                    arrows.append((a6[0], a6[1], math.atan(-dy/dx), right, colors[i]))
                pathRight.append((Path.MOVETO,p2))
                if coordinates[edge[2]]>=breaky:
                    a1=(breakx,(breaky+coordinates[edge[2]])/2)
                    a2=(0, coordinates[edge[2]]-.5)
                    t=approximate(p2, a1, a2, 0.2*width)
                    dx,dy=dydx(p2, a1, a2, t)
                    a3=quad(p2,a1,a2,t)
                    deltaX=breakx-a3[0]
                    deltaY=deltaX*(dy/dx)
                    a4=(breakx, a3[1]+deltaY)
                    pathRight.append((Path.CURVE3,a4))
                    pathRight.append((Path.CURVE3,a3))
                    arrows.append((a3[0], a3[1], math.atan(-dy/dx), right, colors[i]))
                else:
                    if edge[2]>=edge[0]:
                        if coordinates[edge[1]]==breaky:
                            a1=(breakx, breaky+breakx*.75)
                            a2=(breakx/2, breaky+breakx*.75)
                            a3=(breakx/2, breaky+breakx*.25) 
                        else:
                            a1=(breakx, breaky+breakx/2)
                            a2=(breakx/2, breaky+breakx/2)
                            a3=(breakx/2, breaky)
                        a4=(breakx/2, (breaky+coordinates[edge[2]])/2)
                    else:
                        a1=(breakx, breaky+(width-breakx)/2)
                        a2=((width+breakx)/2, breaky+(width-breakx)/2)
                        a3=((width+breakx)/2, breaky)
                        a4=((width+breakx)/2, (breaky+coordinates[edge[2]])/2)
                    a5=(0, coordinates[edge[2]]+.5)
                    pathRight.append((Path.CURVE4, a1))
                    pathRight.append((Path.CURVE4, a2))
                    pathRight.append((Path.CURVE4, a3))
                    t=approximate(a3, a4, a5, 0.2*width)
                    dx,dy=dydx(a3, a4, a5, t)
                    a6=quad(a3,a4,a5,t)
                    deltaX=a3[0]-a6[0]
                    deltaY=deltaX*(dy/dx)
                    a7=(a3[0], a6[1]+deltaY)
                    pathRight.append((Path.CURVE3,a7))
                    pathRight.append((Path.CURVE3,a6))
                    arrows.append((a6[0], a6[1], math.atan(-dy/dx), right, colors[i]))
            codes, verts = zip(*pathRight)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            right.add_patch(patch)
        else:
            p0=(0, coordinates[edge[0]])
            p1=(breakx, (breaky+coordinates[edge[0]])/2)
            p2=(breakx,breaky)
            textPoint=quad(p0,p1,p2,breakx/width)
            left.text(textPoint[0]+.1*width, textPoint[1], cleanProb(edge[3]), fontdict={"color":colors[i]}, horizontalalignment='center', bbox=dict(facecolor=background, alpha=0.7, ec="None"))
            pathLeft=[]
            pathLeft.append((Path.MOVETO, p0))
            pathLeft.append((Path.CURVE3,p1))
            pathLeft.append((Path.CURVE3,p2))
            if edge[1]==edge[2]:
                a1=(breakx,breaky+(coordinates[edge[1]]-breaky)/2)
                a2=(0, coordinates[edge[1]]+.5)
                t=approximate(p2, a1, a2, 0.2*width)
                dx,dy=dydx(p2, a1, a2, t)
                a3=quad(p2,a1,a2,t)
                deltaX=breakx-a3[0]
                deltaY=deltaX*(dy/dx)
                a4=(breakx, a3[1]+deltaY)
                pathLeft.append((Path.CURVE3,a4))
                pathLeft.append((Path.CURVE3,a3))
                arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
                a1=(breakx,coordinates[edge[1]]+.5)
                a2=(0, coordinates[edge[1]]+.5)
                pathLeft.append((Path.MOVETO,p2))
                t=approximate(p2, a1, a2, 0.2*width)
                dx,dy=dydx(p2, a1, a2, t)
                a3=quad(p2,a1,a2,t)
                deltaX=breakx-a3[0]
                deltaY=deltaX*(dy/dx)
                a4=(breakx, a3[1]+deltaY)
                pathLeft.append((Path.CURVE3,a4))
                pathLeft.append((Path.CURVE3,a3))
                arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
            else:
                if coordinates[edge[1]]<=breaky:
                    a1=(breakx,(breaky+coordinates[edge[1]])/2)
                    a2=(0, coordinates[edge[1]]+.5)
                    t=approximate(p2, a1, a2, 0.2*width)
                    dx,dy=dydx(p2, a1, a2, t)
                    a3=quad(p2,a1,a2,t)
                    deltaX=breakx-a3[0]
                    deltaY=deltaX*(dy/dx)
                    a4=(breakx, a3[1]+deltaY)
                    pathLeft.append((Path.CURVE3,a4))
                    pathLeft.append((Path.CURVE3,a3))
                    arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
                else:
                    if edge[1]<=edge[0]:
                        if coordinates[edge[1]]==breaky:
                            a1=(breakx, breaky-breakx*.75)
                            a2=(breakx/2, breaky-breakx*.75)
                            a3=(breakx/2, breaky-breakx*.25) 
                        else:
                            a1=(breakx, breaky-breakx/2)
                            a2=(breakx/2, breaky-breakx/2)
                            a3=(breakx/2, breaky)
                        a4=(breakx/2, (breaky+coordinates[edge[1]])/2)
                    else:
                        a1=(breakx, breaky-(width-breakx)/2)
                        a2=((width+breakx)/2, breaky-(width-breakx)/2)
                        a3=((width+breakx)/2, breaky)
                        a4=((width+breakx)/2, (breaky+coordinates[edge[1]])/2)
                    a5=(0, coordinates[edge[1]]-.5)
                    pathLeft.append((Path.CURVE4, a1))
                    pathLeft.append((Path.CURVE4, a2))
                    pathLeft.append((Path.CURVE4, a3))
                    t=approximate(a3, a4, a5, 0.2*width)
                    dx,dy=dydx(a3, a4, a5, t)
                    a6=quad(a3,a4,a5,t)
                    deltaX=a3[0]-a6[0]
                    deltaY=deltaX*(dy/dx)
                    a7=(a3[0], a6[1]+deltaY)
                    pathLeft.append((Path.CURVE3,a7))
                    pathLeft.append((Path.CURVE3,a6))
                    arrows.append((a6[0], a6[1], math.atan(-dy/dx), left, colors[i]))
                pathLeft.append((Path.MOVETO,p2))
                if coordinates[edge[2]]<=breaky:
                    a1=(breakx,(breaky+coordinates[edge[2]])/2)
                    a2=(0, coordinates[edge[2]]+.5)
                    t=approximate(p2, a1, a2, 0.2*width)
                    dx,dy=dydx(p2, a1, a2, t)
                    a3=quad(p2,a1,a2,t)
                    deltaX=breakx-a3[0]
                    deltaY=deltaX*(dy/dx)
                    a4=(breakx, a3[1]+deltaY)
                    pathLeft.append((Path.CURVE3,a4))
                    pathLeft.append((Path.CURVE3,a3))
                    arrows.append((a3[0], a3[1], math.atan(-dy/dx), left, colors[i]))
                else:
                    if edge[2]<=edge[0]:
                        if coordinates[edge[1]]==breaky:
                            a1=(breakx, breaky-breakx*.75)
                            a2=(breakx/2, breaky-breakx*.75)
                            a3=(breakx/2, breaky-breakx*.25) 
                        else:
                            a1=(breakx, breaky-breakx/2)
                            a2=(breakx/2, breaky-breakx/2)
                            a3=(breakx/2, breaky)
                        a4=(breakx/2, (breaky+coordinates[edge[2]])/2)
                    else:
                        a1=(breakx, breaky-(width-breakx)/2)
                        a2=((width+breakx)/2, breaky-(width-breakx)/2)
                        a3=((width+breakx)/2, breaky)
                        a4=((width+breakx)/2, (breaky+coordinates[edge[2]])/2)
                    a5=(0, coordinates[edge[2]]-.5)
                    pathLeft.append((Path.CURVE4, a1))
                    pathLeft.append((Path.CURVE4, a2))
                    pathLeft.append((Path.CURVE4, a3))
                    t=approximate(a3, a4, a5, 0.2*width)
                    dx,dy=dydx(a3, a4, a5, t)
                    a6=quad(a3,a4,a5,t)
                    deltaX=a3[0]-a6[0]
                    deltaY=deltaX*(dy/dx)
                    a7=(a3[0], a6[1]+deltaY)
                    pathLeft.append((Path.CURVE3,a7))
                    pathLeft.append((Path.CURVE3,a6))
                    arrows.append((a6[0], a6[1], math.atan(-dy/dx), left, colors[i]))
            codes, verts = zip(*pathLeft)
            path = mpath.Path(verts, codes)
            patch = patches.PathPatch(path, linewidth=2, fill=False, aa=True, capstyle="round", edgecolor=colors[i])
            left.add_patch(patch)
    for arrow in arrows:
        makeArrow(*arrow)
    for arrow in arrows:
        makeArrow2(*arrow)
#threshold must be >=1
def generateColors(num, threshold=1):
    #print(num)
    divisions=math.ceil((6*num)**(1/3))
    #print(divisions)
    colors=[]
    frac=1.0/(divisions-1)
    for i in range(divisions):
        for j in range(divisions):
            for k in range(divisions):
                if i*frac+j*frac+k*frac<=threshold:
                    colors.append((i*frac,j*frac,k*frac))
    #print(colors)
    random.shuffle(colors)
    return colors[:num]

def approximate(a,b,c, distance):
    t=.5
    adjustment=0.25
    for i in range(10):
        if (c[0]-((1-t)**2*a[0]+2*(1-t)*t*b[0]+t**2*c[0]))**2+(c[1]-((1-t)**2*a[1]+2*(1-t)*t*b[1]+t**2*c[1]))**2>distance**2:
            t+=adjustment
        else:
            t+=-adjustment
        adjustment=adjustment/2
    return t

def dydx(a,b,c,t):
    return 2*(1-t)*(b[0]-a[0])+2*t*(c[0]-b[0]),2*(1-t)*(b[1]-a[1])+2*t*(c[1]-b[1])

def quad(a,b,c,t):
    return ((1-t)**2*a[0]+2*(1-t)*t*b[0]+t**2*c[0],(1-t)**2*a[1]+2*(1-t)*t*b[1]+t**2*c[1])

def cleanProb(prob):
    s=str(prob)
    if s[0]!="0":
        return str(round(prob,2))
    for i in range(2,len(s)):
        if s[i]!=0:
            return str(round(prob,i+1))
    return s
        

def makeArrow(x, y, angle, axes, color):
    pathArrow=[]
    length=.5
    delta=.4
    pathArrow.append((Path.MOVETO,(x,y)))
    pathArrow.append((Path.LINETO,(x+(length*math.cos(angle+delta)),y-(length*math.sin(angle+delta)))))
    pathArrow.append((Path.LINETO,(x+(length*math.cos(angle-delta)),y-(length*math.sin(angle-delta)))))
    pathArrow.append((Path.LINETO,(x,y)))
    codes, verts = zip(*pathArrow)
    path = mpath.Path(verts, codes)
    patch = patches.PathPatch(path, linewidth=2, fill=True, aa=True, capstyle="round", edgecolor=background, facecolor=color)
    axes.add_patch(patch)
    
def makeArrow2(x, y, angle, axes, color):
    pathArrow=[]
    length=.5
    delta=.4
    pathArrow.append((Path.MOVETO,(x,y)))
    pathArrow.append((Path.LINETO,(x+(length*math.cos(angle+delta)),y-(length*math.sin(angle+delta)))))
    pathArrow.append((Path.LINETO,(x+(length*math.cos(angle-delta)),y-(length*math.sin(angle-delta)))))
    pathArrow.append((Path.LINETO,(x,y)))
    codes, verts = zip(*pathArrow)
    path = mpath.Path(verts, codes)
    patch = patches.PathPatch(path, linewidth=0, fill=True, aa=True, capstyle="round", edgecolor=(1,1,1), facecolor=color)
    axes.add_patch(patch)
    
def makeDiagram():
    # Create the DNAplotlib renderer
    dr = dpl.DNARenderer(linewidth=2.0)
    # Use default renderers and append our custom ones for recombinases
    reg_renderers = dr.std_reg_renderers()
    reg_renderers['Connection'] = flip_arrow
    part_renderers = dr.SBOL_part_renderers()
    part_renderers['RecombinaseSite'] = sbol_recombinase1
    part_renderers['RecombinaseSite2'] = sbol_recombinase2
    part_renderers['RecombinaseSite_round']= sbol_recombinase_round1
    part_renderers['RecombinaseSite_round2']= sbol_recombinase_round2
    part_renderers['RecombinaseSite_arrow']= sbol_recombinase_arrow1
    part_renderers['RecombinaseSite_arrow2']= sbol_recombinase_arrow2
    part_renderers['RecombinaseSite_flying']= sbol_recombinase_flying1
    part_renderers['RecombinaseSite_flying2']= sbol_recombinase_flying2
    part_renderers['RecombinaseSite_block']= sbol_recombinase_block1
    part_renderers['RecombinaseSite_block2']= sbol_recombinase_block2
    part_renderers['RecombinaseSite_nepalese_flag'] = sbol_recombinase_nepalese_flag
    part_renderers['RecombinaseSite_nepalese_flag2'] = sbol_recombinase_nepalese_flag2
    
    
    numLines=36
    width=15
    height=2.5*numLines
    fig = plt.figure(figsize=(width,height))
    gs = gridspec.GridSpec(numLines, 1)
    for i in range(36):
        reader =csv.reader(open("Protein_Math/Circuit"+str(i)+".csv"))
        sequence = [x[0] for x in list(reader)]
        design = grn_dnaplotlib_mapping(sequence)
        ax_dna = plt.subplot(gs[i,0])
        start, end = dr.renderDNA(ax_dna, design, part_renderers, regs=None, reg_renderers=None)
        ax_dna.set_xlim([start, end])
        ax_dna.set_ylim([-30,30])
        ax_dna.set_aspect('equal')
        ax_dna.set_xticks([])
        ax_dna.set_yticks([])
        ax_dna.axis('off')
    fig.savefig("Protein_Math/allCircuits.pdf", transparent=False)
    
makeDiagram()
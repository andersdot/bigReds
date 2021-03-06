#generate pretty figures with proper proportions
import matplotlib as mpl
import matplotlib.pyplot as plt

params = {
    'axes.labelsize' : 9,
    'font.size' : 9,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'text.usetex': False,
    'figure.figsize': [4.5, 4.5]
    }
mpl.rcParams.update(params)

def figsize_and_margins(plotsize,subplots=(1,1),**absolute_margins):
    '''Determine figure size and margins from plot size and absolute margins

       Parameters:
         plotsize: (width, height) of plot area
         subplots: (nrows, ncols) of subplots
         left, right, top, bottom: absolute margins around plot area
         wspace, hspace: width and height spacing between subplots
       Returns:
         size: figure size for figsize argument of figure()
         margins: relative margins dict suitable for subplots_adjust()

       Example: making 2x2 grid of 3" square plots with specific spacings:

       sz, rm = figsize_and_margins((3,3), (2,2), left=1, right=.5,
                                                  top=.5, bottom=1,
                                                  wspace=.5, hspace=.5)
       figure(figsize=sz)
       subplots_adjust(**rm)
       subplot(221); subplot(222)
       subplot(223); subplot(224)
    '''
    #from matplotlib import rcParams
    pw,ph = plotsize
    nr,nc = subplots
    amarg = absolute_margins
    #dictionary for relative margins
    # initialize from rcParams with margins not in amarg
    rmarg = dict((m, mpl.rcParams['figure.subplot.' + m])
                for m in ('left','right','top','bottom','wspace','hspace')
                if m not in amarg
            )
    #subplots_adjust wants wspace and hspace relative to plotsize:
    if 'wspace' in amarg: rmarg['wspace'] = float(amarg['wspace']) / pw
    if 'hspace' in amarg: rmarg['hspace'] = float(amarg['hspace']) / ph
    #in terms of the relative margins:
    #width  * (right - left)
    #    = ncols * plot_width  + (ncols - 1) * wspace * plot_width
    #height * (top - bottom)
    #    = nrows * plot_height + (nrows - 1) * hspace * plot_height
    #solve for width and height, using absolute margins as necessary:

    width  = float((nc + (nc - 1) * rmarg['wspace']) * pw        \
                   + amarg.get('left',0) + amarg.get('right',0)) \
             / (rmarg.get('right',1) - rmarg.get('left',0))

    height = float((nr + (nr - 1) * rmarg['hspace']) * ph        \
                   + amarg.get('top',0) + amarg.get('bottom',0)) \
             / (rmarg.get('top',1) - rmarg.get('bottom',0))

    #now we can get any remaining relative margins
    if 'left'   in amarg: rmarg['left']   =     float(amarg['left'])   / width
    if 'right'  in amarg: rmarg['right']  = 1 - float(amarg['right'])  / width
    if 'top'    in amarg: rmarg['top']    = 1 - float(amarg['top'])    / height
    if 'bottom' in amarg: rmarg['bottom'] =     float(amarg['bottom']) / height
    #return figure size and relative margins
    return (width, height), rmarg

#Example usage: make 2 side-by-side 3" square figures
#from pylab import *
#fsize, margins = figsize_and_margins(plotsize=(3,3),subplots=(1,2))
#figure('My figure', figsize=fsize)
#adjust_subplots(**margins)
#subplot(121)
# ... plot something ...
#subplot(122)
# ... plot something ...
#show()

def makeFigureInstance(x=1, y=1, left=0.5, right=0.25, top=0.25, bottom=0.5, wspace=0.5, hspace=0.5, figureSize=(3,3)):#, figsize=None, fontsize=12):
    sz, rm = figsize_and_margins(figureSize, (y,x), left=left, right=right,
                                           top=top, bottom=bottom,
                                           wspace=wspace, hspace=hspace)
    fig, ax = plt.subplots(y, x, figsize=sz)#, figsize=figsize)
    if (x > 1) or (y > 1): ax = ax.flatten()

    fig.subplots_adjust(**rm)
    #fig.subplots_adjust(left=0.1, right=0.9,
    #                    bottom=0.1, top=0.9,
    #                    wspace=0.4, hspace=0.5)
    #setup_text_plots(fontsize=fontsize, usetex=True)
    return fig, ax

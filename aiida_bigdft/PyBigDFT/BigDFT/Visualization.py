"""
This module has the routines and data structures necessary to allow one to
generate visualizations of a atomic systems. We also include helper routines
for generating colors.

"""


class InlineVisualizer():
    """
    This class allows for a quick viewing of BigDFT systems using the
    py3Dmol package.

    https://pypi.org/project/py3Dmol/

    Attributes:
      xsize (int): the width of the picture in pixels.
      ysize (int): the height of the picture in pixels.
      nrow (int): if present, the number of rows for displaying a grid of
          structures.
      ncol (int): if present, the number of columns for displaying a grid of
          structures.
    """

    def __init__(self, xsize, ysize, nrow=1, ncol=1):
        from py3Dmol import view
        # set linked to False so that each structure can be rotated
        # independently
        self.xyzview = view(width=ncol*xsize, height=nrow*ysize,
                            viewergrid=(nrow, ncol), linked=False)

    def display_system(self, *syslist, **kwargs):
        """
        Display an animation of a sequence of systems. The colordict can be
        used to color each fragment. When only one system is passed it
        will remain still.

        Args:
          syslist (BigDFT.Systems.System): the systems to visualize.
          colordict (dict): a dictionary from fragment ids to hex colors,
             can also be a list of dicts (one for each system) if using a grid.
          field_vals(list): values of the field to decide colors of the keys
          cartoon (bool): set to True to use the cartoon view. This only works
            if atom names and residues are properly defined.
          gridlist (list): if present, defines the row and column indices for
            visualizing multiple systems on a grid.
          show (bool): you can explicitly defer showing.
        """
        from BigDFT.IO import write_pdb
        # py2 workaround
        from sys import version_info
        if version_info[0] < 3:
            from io import BytesIO as StringIO
        else:
            try:
                from io import StringIO
            except ImportError:
                from StringIO import StringIO

        colordict = kwargs.get('colordict')
        field_vals = kwargs.get('field_vals')
        cartoon = kwargs.get('cartoon', False)
        gridlist = kwargs.get('gridlist', None)
        stick_radius = kwargs.get('stick_radius', 0.0)
        stick_color = kwargs.get('stick_color', 'black')
        show = kwargs.get('show', True)

        # Set the default colors
        if colordict is None:
            keyset = []
            for system in syslist:
                keyset += list(system.keys())
            keyset = list(set(keyset))
            colordict = get_colordict(keyset, field_vals=field_vals)

        # Draw each system
        models = ""
        for s, system in enumerate(syslist):
            if type(colordict) is list:
                this_colordict = colordict[s]
            else:
                this_colordict = colordict

            if type(stick_color) is list:
                this_stick_color = stick_color[s]
            else:
                this_stick_color = stick_color

            model = ""
            sval = StringIO()
            write_pdb(system, sval)
            model += "MODEL " + str(s) + "\n"
            model += sval.getvalue()
            model += "ENDMDL\n"

            # in this case we have only one (combined) system to display
            if gridlist is None:
                models += model
            # in this case we just want to display this system at the moment
            else:
                models = model

            # If displaying on a grid, remove all existing models
            if gridlist is not None:
                gx = gridlist[s][0]
                gy = gridlist[s][1]
                viewer = (gx, gy)
                # print(viewer)
                self.xyzview.removeAllModels(viewer=viewer)
            else:
                viewer = (0, 0)

            # in the case of a grid we need to display at each iteration
            # otherwise we only display for the final iteration
            if (gridlist is not None) or (gridlist is None and
                                          s == len(syslist) - 1):
                i = 0
                for fragid, frag in syslist[0].items():
                    if fragid in this_colordict:
                        color = this_colordict[fragid]
                    else:
                        color = 'black'

                    if not cartoon:
                        if stick_radius > 0.0:
                            self.xyzview.addModelsAsFrames(models, "pdb",
                                                           {"keepH": "true"},
                                                           viewer=viewer)
                            self.xyzview.setStyle({'model': -1},
                                                  {"stick":
                                                  {'radius': stick_radius,
                                                   'color': this_stick_color}},
                                                  viewer=viewer)
                        else:
                            self.xyzview.addModelsAsFrames(models, "pdb",
                                                           {"keepH": "true"},
                                                           viewer=viewer)
                            self.xyzview.setStyle({'model': -1},
                                                  {"line": {'color': 'black'}},
                                                  viewer=viewer)

                    for at in frag:
                        if cartoon:
                            self.xyzview.addModelsAsFrames(models,
                                                           viewer=viewer)
                            self.xyzview.setStyle({'model': -1,
                                                   'serial': i + 1},
                                                  {"cartoon":
                                                  {'color': color}},
                                                  viewer=viewer)
                        else:

                            self.xyzview.addModelsAsFrames(models, "pdb",
                                                           {"keepH": "true"},
                                                           viewer=viewer)
                            self.xyzview.setStyle({'model': -1, 'serial': i+1},
                                                  {"sphere": {'scale': 0.2,
                                                              'color': color}},
                                                  viewer=viewer)
                        i += 1

            if gridlist is not None:
                self.xyzview.render()
        self.display_cell(syslist[0].cell)

        # Finish Drawing
        if len(syslist) > 1 and gridlist is None:
            self.xyzview.animate({'loop': "forward", 'interval': 1000})
        self.xyzview.zoomTo()

        if show:
            self.xyzview.show()

    def display_cell(self, cell):
        if cell is None:
            return
        self.xyzview.addUnitCell({"model": -1},
                                 {"box": {"color": "black"},
                                 "alabel": "", "blabel": "", "clabel": ""})


class VSimVisualizer():
    def __init__(self, filename, xsize=600, ysize=600):
        from gi.repository import v_sim, GLib, Gtk
        self.win = v_sim.UiRenderingWindow.new(xsize, ysize, True, True)
        self.loop = GLib.MainLoop.new(None, False)
        self.win.connect_object('destroy', GLib.MainLoop.quit, self.loop)
        self.main = Gtk.Window.new(Gtk.WindowType.TOPLEVEL)
        self.main.add(self.win)
        self.main.show_all()
        self._set_file(filename)

    def _set_file(self, filename):
        from gi.repository import v_sim
        self.data = v_sim.DataAtomic.new(filename, None)
        self.data.load(0, None)
        self.win.getGlScene().setData(self.data)
        # v_sim.basic_parseConfigFiles()

    def show(self):
        # from gi.repository import GLib
        self.loop.run()

    def colorize_by_fragments(self):
        from gi.repository import v_sim
        self._scene_colorizer(v_sim, self.data, self.win.getGlScene())

    def _scene_colorizer(self, v_sim, data, scene):
        nodes = scene.getNodes()
        frag = data.getNodeProperties("Fragment")
        c = v_sim.DataColorizerFragment.new()
        c.setNodeModel(frag)
        nodes.pushColorizer(c)
        c.setActive(True)

    def colorizer_script(self, filename):
        import inspect
        towrite = ['scene = ' +
                   'v_sim.UiMainClass.getDefaultRendering().getGlScene()',
                   'data = scene.getData()']
        for line in inspect.getsource(self._scene_colorizer).split('\n')[1:]:
            towrite.append(line.lstrip(' '))
        f = open(filename, 'w')
        for line in towrite:
            f.write(line + '\n')
        f.close()


class VMDGenerator():
    """
    This class contains the routines you would use for visualization of
    a system using the VMD program.

    Attributes:
      representation (str): the vmd representation to draw with.
        https://www.ks.uiuc.edu/Research/vmd/allversions/repimages/#representations
      color (int): the default color to draw with.
    """

    def __init__(self, representation="CPK", color=16):
        self.representation = representation
        self.color = color

    def visualize_fragments(self, system, scriptfile, geomfile,
                            fragcolors=None):
        """
        This generates a script for visualizing the fragmentation of a
        system using VMD.

        Args:
          system (BigDFT.Systems.System): the system to visualize.
          scriptfile (str): the name of the file to write the vmd script
            to (usually has extension .tcl)
          geomfile (str): the filename for where to write an xyz file
            of the system.
          fragcolors (dict): optionally, a dictionary from fragment ids to
            fragment colors. Colors are integers between 0 and 32.
        """
        from BigDFT.Fragments import Fragment
        from BigDFT.IO import XYZWriter

        # To create the XYZ file, we first make one big fragment.
        geomorder = Fragment()
        for fragid, frag in system.items():
            geomorder += frag

        # Then write it to file.
        with XYZWriter(geomfile, len(geomorder)) as ofile:
            for at in geomorder:
                ofile.write(at)

        # Get the matching so we can write the correct atom indices.
        matching = system.compute_matching(geomorder)

        # If fragcolors is not specified, we will generate it ourselves.
        if fragcolors is None:
            fragcolors = {}
            for i, s in enumerate(system):
                # 16 is black, which we will reserve.
                if i % 32 == 16:
                    c = str(32)
                else:
                    c = str(i % 32)
                fragcolors[s] = c

        # The header of the script file draws the whole system in black.
        outstr = self._get_default_header(geomfile)

        # This part colors individual fragments.
        modid = 1
        for fragid, frag in system.items():
            if fragid not in fragcolors:
                continue
            outstr += "mol addrep 0\n"
            outstr += """mol modselect """ + str(modid) + """ 0 index """
            outstr += " ".join([str(x) for x in matching[fragid]])
            outstr += "\n"
            outstr += """mol modcolor """
            outstr += str(modid) + """ 0 ColorID """ + \
                str(fragcolors[fragid]) + """\n"""
            modid += 1

        # Finally, write to file.
        with open(scriptfile, "w") as ofile:
            ofile.write(outstr)

    def _get_default_header(self, geomfile):
        outstr = """mol default style """+self.representation+"""\n"""
        outstr += """mol new """ + geomfile + "\n"
        outstr += """mol modcolor 0 0 ColorID """+str(self.color)+"""\n"""

        return outstr


def get_distinct_colors(keys, name="tab20", fuzz=True):
    """
    This generates a dictionary of distinct colors based on a matplotlib
    colormap.

    Args:
        keys (list): a list of keys.
        name (str): the name of the matplotlib colormap to use.
        fuzz (bool): some color maps (included tab20) only have a distinct
            set of colors. The fuzz option adds increased randomness to make
            up for this.

    Returns:
        (dict): a dictionary mapping matplotlib keys to RGB colors.
    """
    from matplotlib import pyplot as plt
    from numpy import linspace
    from random import sample, uniform

    def fuzz(x):
        fuzzed = (x[i] + uniform(-0.1, 0.1) for i in range(3))
        return [max(0.0, min(y, 1.0)) for y in fuzzed]

    # Map using a color map
    cmap = plt.get_cmap(name)
    xvals = linspace(0, 1, len(keys))
    vals = {x: list(cmap(xvals[i])) for i, x in
            enumerate(sample(keys, len(keys)))}

    # Add some random fuzz
    if fuzz:
        vals = {x: fuzz(y) for x, y in vals.items()}

    return vals


def get_colordict(keys, field_vals=None, colorcode='html'):
    """
    Build a dictionary of colors for each of the keys.
    If the field_dict is provided, order the colors of the keys in terms
    of the sorting of the filed values

    Args:
        keys (list): keys of the color dictionary
        field_vals(list) : values of the field to decide the colors of the keys
        colorcode (str): the string of the colorcode.
             Can be 'html' or 'rgb' if field_dict is absent,
             otherwise it represent the colormap of matplotlib.
             Default is 'seismic' for diverging values
             (field_dict has negative data), otherwise 'Reds'
    Returns:
        dict: the dictionary of the keys, and the corresponding colors
    """
    import numpy as np
    from matplotlib.pyplot import get_cmap
    cvals = _find_colours(len(keys))
    if field_vals is None:
        colorkey = 'indices' if colorcode == 'rgb' else colorcode
        cvals = cvals[colorkey]
    else:
        compressed_values = np.array(
            [val for val in field_vals if not np.isnan(val)])
        center_to_zero = any(compressed_values < 0.0)
        mx = np.max(compressed_values)
        mn = np.min(compressed_values)
        if center_to_zero:
            shift = min(mn, -mx)
            top = max(mx, -mn)
        else:
            shift = mn
            top = mx
        compressed_values -= shift
        if top > shift:
            compressed_values /= (top-shift)
        if colorcode == 'html':
            colorkey = 'seismic' if center_to_zero else 'Reds'
        else:
            colorkey = colorcode
        cmap = get_cmap(colorkey)
        cvals = []
        icv = 0
        for v in field_vals:
            if np.isnan(v):
                cvals.append('None')
            else:
                y = compressed_values[icv]
                cvals.append(_rgb_to_html(
                    tuple(map(int, np.array(cmap(y)[0:3])*255))))
                icv += 1
        # cvals = [_rgb_to_html(tuple(map(int, np.array(cmap(y)[0:3])*255)))
        #          for y in compressed_values]

    return {x: y for x, y in zip(keys, cvals)}


def _find_colours(n):
    """
    A helper routine for finding a list of N colors in varous formats.

    Args:
      n (int): the number of colors to generate.
           Cycles colours after 1408 values

    Returns:
      (dict): a dictionary describing the colors. The key `html` contains
      the hex values values, and `indices` contains a list of RGB values.
    """
    colours = {}
    colours['html'] = list()
    colours['indices'] = list()
    ind3d = [0] * 3
    for i in range(n):
        ind1d = 128+i*int(1408/n)
        if (ind1d < 256):
            ind3d[0] = ind1d
            ind3d[1] = 0
            ind3d[2] = 0
        elif (ind1d < 512):
            ind3d[0] = 255
            ind3d[1] = ind1d-256
            ind3d[2] = 0
        elif (ind1d < 768):
            ind3d[0] = 768-ind1d-1
            ind3d[1] = 255
            ind3d[2] = 0
        elif (ind1d < 1024):
            ind3d[0] = 0
            ind3d[1] = 255
            ind3d[2] = ind1d-768
        elif (ind1d < 1280):
            ind3d[0] = 0
            ind3d[1] = 1280-ind1d-1
            ind3d[2] = 255
        elif (ind1d <= 1536):
            ind3d[0] = ind1d-1280
            ind3d[1] = 0
            ind3d[2] = 255

        colours['indices'].append(ind3d[:])
        html = _rgb_to_html(ind3d)

        colours['html'].append(html)

    return colours


def _rgb_to_html(rgb):
    html = [0] * 6
    htmls = [0] * 6
    html[0] = int(rgb[0]/16)
    html[1] = (rgb[0] % 16)
    html[2] = int(rgb[1]/16)
    html[3] = (rgb[1] % 16)
    html[4] = int(rgb[2]/16)
    html[5] = (rgb[2] % 16)
    for j in range(6):
        if (html[j] == 10):
            htmls[j] = 'A'
        elif (html[j] == 11):
            htmls[j] = 'B'
        elif (html[j] == 12):
            htmls[j] = 'C'
        elif (html[j] == 13):
            htmls[j] = 'D'
        elif (html[j] == 14):
            htmls[j] = 'E'
        elif (html[j] == 15):
            htmls[j] = 'F'
        else:
            htmls[j] = html[j]
        html_string = '#'+str(htmls[0])+str(htmls[1]) + \
            str(htmls[2])+str(htmls[3])+str(htmls[4])+str(htmls[5])
    return html_string


def _example():
    """Visualization Example"""
    from BigDFT.Systems import System
    from BigDFT.Fragments import Fragment
    from BigDFT.IO import XYZReader

    # Read in a system.
    sys = System()
    with XYZReader("SiO") as ifile:
        for i, at in enumerate(ifile):
            sys["FRA:"+str(i)] = Fragment([at])

    # Display the system.
    viz = InlineVisualizer(400, 300)
    viz.display_system(sys)

    # Change the colors
    colordict = get_distinct_colors(list(sys))
    viz.display_system(sys, colordict=colordict)


if __name__ == "__main__":
    _example()

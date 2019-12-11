"""
This module has the routines and data structures necessary to allow one to
generate visualizations of a atomic systems.

"""


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
          system (Fragments.System): the system to visualize.
          scriptfile (str): the name of the file to write the vmd script
            to (usually has extension .tcl)
          geomfile (str): the filename for where to write an xyz file
            of the system.
          fragcolors (dict): optionally, a dictionary from fragment ids to
            fragment colors. Colors are integers between 0 and 32.
        """
        from BigDFT.Fragments import Fragment
        from BigDFT.XYZ import XYZWriter

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

    def visualize_qmmm(self, system, subsystem, target, scriptfile,
                       system_geomfile, subsystem_geomfile):
        """
        This routine generates a VMD script for visualizing a QM/MM setup.

        Args:
          system (Fragments.System): the full system to visualize.
          subsystem (Fragments.System): the qm part to visualize.
          target (str): the id of the target fragment of the QM/MM calculation.
          scriptfile (str): the name of the file to write the vmd script
            to (usually has extension .tcl)
          system_geomfile (str): the filename for where to write an xyz file
            of the system.
          subsystem_geomfile (str): the filename for where to write an xyz file
            of the subsystem.
        """
        from BigDFT.Fragments import Fragment
        from BigDFT.XYZ import XYZWriter

        # To create the XYZ file, we first make one big fragment.
        systemorder = Fragment()
        for fragid, frag in system.items():
            systemorder += frag
        subsystemorder = Fragment()
        for fragid, frag in subsystem.items():
            subsystemorder += frag

        # Then write it to file.
        with XYZWriter(system_geomfile, len(systemorder)) as ofile:
            for at in systemorder:
                ofile.write(at)
        with XYZWriter(subsystem_geomfile, len(subsystemorder)) as ofile:
            for at in subsystemorder:
                ofile.write(at)

        # The header of the script file draws the whole system in black.
        outstr = self._get_default_header(system_geomfile)

        # Draw the QM portion of the system
        matching = system.compute_matching(systemorder)
        modid = 1
        for fragid, frag in subsystem.items():
            if fragid == target:
                color = 4
            else:
                color = 0
            outstr += "mol addrep 0\n"
            outstr += """mol modselect """ + str(modid) + """ 0 index """
            outstr += " ".join([str(x) for x in matching[fragid]])
            outstr += "\n"
            outstr += """mol modcolor """
            outstr += str(modid) + """ 0 ColorID """ + str(color) + """\n"""
            modid += 1

        # Finally, write to file.
        with open(scriptfile, "w") as ofile:
            ofile.write(outstr)

    def _get_default_header(self, geomfile):
        outstr = """mol default style """+self.representation+"""\n"""
        outstr += """mol new """ + geomfile + "\n"
        outstr += """mol modcolor 0 0 ColorID """+str(self.color)+"""\n"""

        return outstr

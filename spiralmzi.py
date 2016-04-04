from ipkiss3 import all as i3
TECH = i3.TECH
from ipcore.caching.cache import cache
from picazzo3.traces.wire_wg import WireWaveguideTemplate
from picazzo3.wg.dircoup import BendDirectionalCoupler
from numpy import all as np
from picazzo3.wg.spirals import DoubleSpiralWithInCouplingRounded 
#from picazzo3.wg.spirals import DoubleSpiralRounded 
from reme import *
from reme_ipkiss.pysimul.reme_view import RemeView

#a class for directional coupler, caphemodel is based on reme results of coupler calculation
class BDCoupler(i3.PCell):
    _name_prefix = "COUPLER"
    simdata = i3.DefinitionProperty(default=0,doc="reme simulation data")
    class Netlist(i3.NetlistView):
        def _generate_terms(self, terms):
            terms += i3.OpticalTerm(name="in1")
            terms += i3.OpticalTerm(name="in2")
            terms += i3.OpticalTerm(name="out1")
            terms += i3.OpticalTerm(name="out2")
            return terms

    class CapheModel(i3.CapheSModelView):
        def _calculate_S(self, environment, term1, term2, mode1, mode2): 
            straight = (self.simdata.s0000_f(environment.wavelength))
            cross = (self.simdata.s1000_f(environment.wavelength))

            if ((term1.name =='in1' and term2.name=='out1') or
               (term1.name=='in2' and term2.name =='out2') or
               (term1.name=='out1' and term2.name=='in1') or
               (term1.name=='out2' and term2.name=='in2')):
                return straight

            if ((term1.name=='in1' and term2.name=='out2') or
               (term1.name=='in2' and term2.name=='out1') or
               (term1.name=='out1' and term2.name=='in2') or
               (term1.name=='out2' and term2.name=='in1')):
                return cross
            return 0

#a class for MZI arms
class SpiralArm(i3.PCell):
    _name_prefix = "SPIRAL_ARM"
    simdata = i3.DefinitionProperty(default=0,doc="reme simulation data")
    loss_dB_cm = i3.NumberProperty(default=0, doc="waveguide loss")
    spiral_length = i3.PositiveNumberProperty(required=True, doc="length of spiral arm part")
    left_length = i3.PositiveNumberProperty(required=True, doc="length of the left of connecting arm part")
    right_length = i3.PositiveNumberProperty(required=True, doc="length of the right of connecting arm part")

    class Netlist(i3.NetlistView):
        def _generate_terms(self, terms):
            terms += i3.OpticalTerm(name="in")
            terms += i3.OpticalTerm(name="out")
            return terms

    class CapheModel(i3.CapheSModelView):
        def _calculate_S(self, environment, term1, term2, mode1, mode2):
            arm_length = self.spiral_length + self.left_length + self.right_length #total length of one arm
            phase = (2*np.pi/environment.wavelength) * self.simdata.fstraight(environment.wavelength) * (arm_length)
            loss_dB_m = self.loss_dB_cm * 100
            loss = 10 ** ((-loss_dB_m * arm_length * 1e-6)/20.0)
            if((term1.name=='in' and term2.name=='out') or (term1.name=='out' and term2.name=='in')):
                return np.exp(1j * phase) * loss
            return 0
        
# class for MZI
class SpiralMZI(i3.PCell):
    _name_prefix = "SPIRAL_MZI"
    wg_template = i3.WaveguideTemplateProperty(default=i3.TECH.PCELLS.WG.DEFAULT)
    wg_dc_template = i3.WaveguideTemplateProperty(default=i3.TECH.PCELLS.WG.DEFAULT)
    wgs = i3.ChildCellListProperty()
    loss_dB_cm = i3.NumberProperty(default=0, doc="waveguide loss")
    splitter = i3.ChildCellProperty()
    combiner = i3.ChildCellProperty()
    spiral_north = i3.ChildCellProperty()
    spiral_south = i3.ChildCellProperty()
    dc_coupler_length = i3.PositiveNumberProperty(required=True,doc="coupler length of directional coupler") 
    dc_gap = i3.PositiveNumberProperty(required=True,doc="gap of directional coupler")  
    radius = i3.PositiveNumberProperty(required=True,doc="bend radius of waveguides") 
    # 2. We define a waveguide template and the waveguide cells 
    def _default_wgs(self):
        wg1 = i3.RoundedWaveguide(name="wg1", trace_template=self.wg_template)
        wg2 = i3.RoundedWaveguide(name="wg2", trace_template=self.wg_template)
        wg3 = i3.RoundedWaveguide(name="wg3", trace_template=self.wg_template)
        wg4 = i3.RoundedWaveguide(name="wg4", trace_template=self.wg_template)
        return wg1, wg2, wg3, wg4
    
    # I take "direction coupler" in Ipkiss library as default splitter and combiner just to exploit the layout view in the library 
    def _default_splitter(self):
        return BendDirectionalCoupler(coupler_length=self.dc_coupler_length)

    def _default_combiner(self):
        return BendDirectionalCoupler(coupler_length=self.dc_coupler_length)
    
    # For the same reason, I take "spiral" in Ipkiss library to make use of its default layout view
    def _default_spiral_north(self):
        return DoubleSpiralWithInCouplingRounded(n_o_loops=3)

    def _default_spiral_south(self):
        return DoubleSpiralWithInCouplingRounded(n_o_loops=2)

    class RemeSimulation(RemeView):
        # As you said, I use rounded waveguide with bend radius > 10um so that I can use solver for straigh guide results

        start = i3.PositiveNumberProperty(default=1.51,doc="start wavelength")
        stop = i3.PositiveNumberProperty(default=1.6,doc="stop wavelength")
        coupler_sampling_x = i3.PositiveNumberProperty(default=501,doc="coupler sampling x")
        coupler_sampling_y = i3.PositiveNumberProperty(default=501,doc="coupler sampling y")
        coupler_step = i3.PositiveNumberProperty(default=51,doc="coupler step")
        coupler_overhang = i3.PositiveNumberProperty(default=10,doc="addition coupler length")

        straight_part_dc = i3.LockedProperty()
        def _default_straight_part_dc(self):
            rwg = self.cell.wg_dc_template.views["rememodemodel"].rwg
            straight_part_dc = FMMStraight(rwg)
            straight_part_dc.set_top_boundary(PMC)
            straight_part_dc.set_bottom_boundary(PMC)
            straight_part_dc.set_left_boundary(PEC)
            straight_part_dc.set_right_boundary(PEC)
            return straight_part_dc

        bent_part_dc = i3.LockedProperty()
        def _default_bent_part_dc(self):
            rwg = self.cell.wg_dc_template.views["rememodemodel"].rwg
            bent_part_dc = FMMBent(rwg)
            bent_part_dc.set_radius(self.radius*1e-6)
            bent_part_dc.set_top_boundary(PMC)
            bent_part_dc.set_bottom_boundary(PMC)
            return bent_part_dc

        coupler_dc = i3.LockedProperty()
        def _default_coupler_dc(self):
            coupler_dc = Coupler(2)
            coupler_dc.add_straight_guide(0,self.straight_part_dc)
            coupler_dc.add_bent_guide(0,self.bent_part_dc)
            coupler_dc.add_straight_guide(1,self.straight_part_dc,self.cell.wg_dc_template.views['layout'].core_width * 1e-6 + 
                                         self.dc_gap * 1e-6)
            coupler_dc.add_bent_guide(1,self.bent_part_dc,self.cell.wg_dc_template.views['layout'].core_width * 1e-6 + 
                                         self.dc_gap * 1e-6)
            coupler_dc.set_coupler_length((self.dc_coupler_length+self.coupler_overhang) * 1e-6)
            coupler_dc.set_straight_length(0,self.dc_coupler_length * 1e-6)
            coupler_dc.set_sampling(int(self.coupler_sampling_x),int(self.coupler_sampling_y))
            coupler_dc.set_num_steps(int(self.coupler_step))
            return coupler_dc

        straight_guide = i3.LockedProperty()
        def _default_straight_guide(self):
            rwg = self.cell.wg_template.views["rememodemodel"].rwg
            straight_guide = FMMStraight(rwg)
            straight_guide.set_top_boundary(PMC)
            straight_guide.set_bottom_boundary(PMC)
            straight_guide.set_left_boundary(PEC)
            straight_guide.set_right_boundary(PEC)
            return straight_guide
        
        def coupler_dc_calculation(self, wavelength):
           set_wavelength(wavelength*1e-6)
           self.coupler_dc.calculate()
           svalues = dict(s0000=self.coupler_dc.get_S_reduced(0,0,0,0),
                          s0010=self.coupler_dc.get_S_reduced(0,0,1,0),
                          s1010=self.coupler_dc.get_S_reduced(1,0,1,0),
                          s1000=self.coupler_dc.get_S_reduced(1,0,0,0))
           return svalues

        def other_calculations(self, wavelength):
           set_wavelength(wavelength*1e-6)
           self.straight_guide.clear_mode_list()    
           self.straight_part_dc.clear_mode_list()    
           self.bent_part_dc.clear_mode_list()

           self.straight_guide.find_modes(1)
           self.straight_part_dc.find_modes(1)
           self.bent_part_dc.find_mode((self.straight_part_dc.get_mode_effective_index(0) * 2.0*np.pi/get_wavelength() * (self.radius*1e-6)))
           return 0

        def calculation(self): 
           self.mid =  self.start + (self.stop-self.start)*0.5

           self.other_calculations(self.start)
           #start_n_bent = (self.bent_guide.get_mode_effective_index(0)/
           #                  (2.0*np.pi/get_wavelength()*(self.radius*1e-6)))
           start_n_straight = self.straight_guide.get_mode_effective_index(0)
           self.start_cp = self.coupler_dc_calculation(self.start)
  
           self.other_calculations(self.mid)
           #mid_n_bent = (self.bent_guide.get_mode_effective_index(0)/
           #                  (2.0*np.pi/get_wavelength()*(self.radius*1e-6)))
           mid_n_straight = self.straight_guide.get_mode_effective_index(0)
           self.mid_cp = self.coupler_dc_calculation(self.mid)

           self.other_calculations(self.stop)
           #stop_n_bent = (self.bent_guide.get_mode_effective_index(0)/
           #                  (2.0*np.pi/get_wavelength()*(self.radius*1e-6)))
           stop_n_straight = self.straight_guide.get_mode_effective_index(0)
           self.stop_cp = self.coupler_dc_calculation(self.stop)

           #polynomial curve fitting
           #z = np.polyfit([self.start, self.mid, self.stop], [start_n_bent, mid_n_bent, stop_n_bent], 2)
           #self.fbent = np.poly1d(z)
           z = np.polyfit([self.start, self.mid, self.stop], [start_n_straight, mid_n_straight, stop_n_straight], 2)
           self.fstraight = np.poly1d(z)

           z = np.polyfit([self.start, self.mid, self.stop], [self.start_cp['s0000'], self.mid_cp['s0000'], self.stop_cp['s0000']], 2)
           self.s0000_f = np.poly1d(z)

           z = np.polyfit([self.start, self.mid, self.stop], [self.start_cp['s1000'], self.mid_cp['s1000'], self.stop_cp['s1000']], 2)
           self.s1000_f = np.poly1d(z)

           #z = np.polyfit([self.start, self.mid, self.stop], [self.start_cp['s1010'], self.mid_cp['s1010'], self.stop_cp['s1010']], 2)
           #self.s1010_f = np.poly1d(z)

           #z = np.polyfit([self.start, self.mid, self.stop], [self.start_cp['s0010'], self.mid_cp['s0010'], self.stop_cp['s0010']], 2)
           #self.s0010_f = np.poly1d(z)
           return 0

    class Layout(i3.LayoutView):
        def _default_wgs(self):
            from ipkiss.plugins.photonics.routing.manhattan import RouteManhattan
            splitter, combiner, spiral_north, spiral_south = self._get_components()
            wg1_cell, wg2_cell, wg3_cell, wg4_cell = self.cell.wgs

            wg1_layout = wg1_cell.get_default_view(i3.LayoutView)
            wg1_layout.set(shape=RouteManhattan(input_port=splitter.ports["out2"], output_port=spiral_north.ports["in"]),
                          bend_radius=self.radius
                          )
            wg2_layout = wg2_cell.get_default_view(i3.LayoutView)
            wg2_layout.set(shape=RouteManhattan(input_port=spiral_north.ports["out"], output_port=combiner.ports["in2"]),
                          bend_radius=self.radius
                          )
            wg3_layout = wg3_cell.get_default_view(i3.LayoutView)
            wg3_layout.set(shape=RouteManhattan(input_port=splitter.ports["out1"], output_port=spiral_south.ports["in"]), 
                          bend_radius=self.radius
                          )
            wg4_layout = wg4_cell.get_default_view(i3.LayoutView)
            wg4_layout.set(shape=RouteManhattan(input_port=spiral_south.ports["out"], output_port=combiner.ports["in1"]), 
                          bend_radius=self.radius
                          )
            return wg1_layout, wg2_layout, wg3_layout, wg4_layout
        
        def get_components_layout_view(self):
            core_width = self.wg_template.core_width
            #print "core width",core_width
            splitter_layout = self.cell.splitter.Layout(coupler_spacing=core_width + self.dc_gap, bend_radius=self.radius)
            combiner_layout = self.cell.splitter.Layout(coupler_spacing=core_width + self.dc_gap, bend_radius=self.radius)
            spiral_south_layout = self.cell.spiral_south.Layout(bend_radius=self.radius)
            spiral_north_layout = self.cell.spiral_north.Layout(bend_radius=self.radius)
            return splitter_layout, combiner_layout, spiral_south_layout, spiral_north_layout 

        def _get_components(self):
            splitter_layout, combiner_layout, spiral_south_layout, spiral_north_layout = self.get_components_layout_view() 
            splitter_transformation = i3.Translation((0.0, 0.0))
            combiner_transformation = i3.Translation((260.0,0.0))
            spiral_north_transformation = i3.Translation((80.0,80.0)) 
            spiral_south_transformation = i3.Translation((80.0,-80.0)) + i3.VMirror(-80.0)
            splitter = i3.SRef(reference=splitter_layout, name="splitter", transformation=splitter_transformation)
            combiner = i3.SRef(reference=combiner_layout, name="combiner", transformation=combiner_transformation)
            spiral_south = i3.SRef(reference=spiral_south_layout, name="spiral_south", transformation=spiral_south_transformation)
            spiral_north = i3.SRef(reference=spiral_north_layout, name="spiral_north", transformation=spiral_north_transformation)
            return splitter, combiner, spiral_north, spiral_south

        def _generate_elements(self, elems):
            wg1_layout, wg2_layout, wg3_layout, wg4_layout = self.wgs
            elems += self._get_components()
            elems += i3.SRef(reference=wg1_layout, name="wg1") 
            elems += i3.SRef(reference=wg2_layout, name="wg2") 
            elems += i3.SRef(reference=wg3_layout, name="wg3")
            elems += i3.SRef(reference=wg4_layout, name="wg4") 
            return elems

        def _generate_ports(self, prts):
            prts += self.splitter.ports["in1"].modified_copy(name="in_1")
            prts += self.splitter.ports["in2"].modified_copy(name="in_2")
            prts += self.combiner.ports["out1"].modified_copy(name="out_1")
            prts += self.combiner.ports["out2"].modified_copy(name="out_2")
            return prts

        def get_arm1_length(self):
            wg1_layout, wg2_layout, wg3_layout, wg4_layout = self.wgs
            return wg1_layout.trace_length() 

        def get_arm2_length(self):
            wg1_layout, wg2_layout, wg3_layout, wg4_layout = self.wgs
            return wg2_layout.trace_length() 

        def get_arm3_length(self):
            wg1_layout, wg2_layout, wg3_layout, wg4_layout = self.wgs
            return wg3_layout.trace_length() 

        def get_arm4_length(self):
            wg1_layout, wg2_layout, wg3_layout, wg4_layout = self.wgs
            return wg4_layout.trace_length() 

        def get_spiral_north_length(self):
            splitter_layout, combiner_layout, spiral_south_layout, spiral_north_layout = self.get_components_layout_view() 
            return spiral_north_layout.trace_length() 

        def get_spiral_south_length(self):
            splitter_layout, combiner_layout, spiral_south_layout, spiral_north_layout = self.get_components_layout_view() 
            return spiral_south_layout.trace_length() 

    class Netlist(i3.NetlistView):
        def _generate_terms(self, terms):
            terms += i3.OpticalTerm(name="in_1")
            terms += i3.OpticalTerm(name="in_2")
            terms += i3.OpticalTerm(name="out_1")
            terms += i3.OpticalTerm(name="out_2")
            return terms

        def _generate_instances(self, insts):
            splitter_mzi = BDCoupler(simdata=self.cell.views['remesimulation'])
            combiner_mzi = splitter_mzi 
            north_arm_mzi = SpiralArm(loss_dB_cm=self.loss_dB_cm, simdata=self.cell.views['remesimulation'],
                                      spiral_length=self.cell.views['layout'].get_spiral_north_length(),
                                      left_length=self.cell.views['layout'].get_arm1_length(),
                                      right_length=self.cell.views['layout'].get_arm2_length())

            south_arm_mzi = SpiralArm(loss_dB_cm=self.loss_dB_cm, simdata=self.cell.views['remesimulation'],
                                      spiral_length=self.cell.views['layout'].get_spiral_south_length(),
                                      left_length=self.cell.views['layout'].get_arm3_length(),
                                      right_length=self.cell.views['layout'].get_arm4_length())

            #print "spiral north length",self.cell.views['layout'].get_spiral_north_length()
            #print "spiral south length",self.cell.views['layout'].get_spiral_south_length()
            #print "arm1 length",self.cell.views['layout'].get_arm1_length()
            #print "arm2 length",self.cell.views['layout'].get_arm2_length()
            #print "arm3 length",self.cell.views['layout'].get_arm3_length()
            #print "arm4 length",self.cell.views['layout'].get_arm4_length()
            print "north arm total length",self.cell.views['layout'].get_spiral_north_length()+self.cell.views['layout'].get_arm1_length()+self.cell.views['layout'].get_arm1_length()

            print "south arm total length",self.cell.views['layout'].get_spiral_south_length()+self.cell.views['layout'].get_arm3_length()+self.cell.views['layout'].get_arm4_length()

            insts += i3.Instance(name="splitter",reference=splitter_mzi)
            insts += i3.Instance(name="combiner",reference=combiner_mzi)
            insts += i3.Instance(name="north_arm",reference=north_arm_mzi)
            insts += i3.Instance(name="south_arm",reference=south_arm_mzi)
            return insts

        def _generate_nets(self, nets):
            nets += i3.OpticalLink(name="li11", term1=self.terms["in_1"], term2=self.instances["splitter"].terms["in1"])
            nets += i3.OpticalLink(name="li22", term1=self.terms["in_2"], term2=self.instances["splitter"].terms["in2"])
            nets += i3.OpticalLink(name="lo11", term1=self.instances["combiner"].terms["out1"], term2=self.terms["out_1"]) 
            nets += i3.OpticalLink(name="lo22", term1=self.instances["combiner"].terms["out2"], term2=self.terms["out_2"]) 

            #nets += i3.OpticalLink(name="li11", term1=self.terms["in_1"], term2=self.instances["splitter"].terms["in1"])
            #nets += i3.OpticalLink(name="li22", term1=self.terms["in_2"], term2=self.instances["splitter"].terms["in2"])
            #nets += i3.OpticalLink(name="lo11", term1=self.instances["splitter"].terms["out1"], term2=self.terms["out_1"]) 
            #nets += i3.OpticalLink(name="lo22", term1=self.instances["splitter"].terms["out2"], term2=self.terms["out_2"]) 

            #nets += i3.OpticalLink(name="l1", term1=self.instances["splitter"].terms["out2"], term2=self.instances["combiner"].terms["in2"]) 
            #nets += i3.OpticalLink(name="l2", term1=self.instances["splitter"].terms["out1"], term2=self.instances["combiner"].terms["in1"]) 

            nets += i3.OpticalLink(name="l1", term1=self.instances["splitter"].terms["out2"], term2=self.instances["south_arm"].terms["in"]) 
            nets += i3.OpticalLink(name="l2", term1=self.instances["splitter"].terms["out1"], term2=self.instances["north_arm"].terms["in"]) 
            nets += i3.OpticalLink(name="l3", term1=self.instances["combiner"].terms["in2"], term2=self.instances["south_arm"].terms["out"]) 
            nets += i3.OpticalLink(name="l4", term1=self.instances["combiner"].terms["in1"], term2=self.instances["north_arm"].terms["out"]) 
            return nets

    class CapheModel(i3.CapheModelFromNetlist):
        pass



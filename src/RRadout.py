### Revit to Radiance scene file export module
### 2016 Georg Mischler
###
### Version 0.2 - prototype
###
###############################################################################
#
# The MIT License (MIT)
#
# Copyright (c) 2016 Georg Mischler, Munich, Germany
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
###############################################################################

__software__ = 'RRadout'
__version__ = '0.2'

import sys
import nt as os
import time
import clr
import System
clr.AddReferenceByName("RevitAPI.dll");
clr.AddReferenceByName("RevitAPIUI.dll");

from Autodesk.Revit import *
from Autodesk.Revit.UI import *
from Autodesk.Revit.UI.Macros import *
from Autodesk.Revit.DB import *
from Autodesk.Revit.DB.Architecture import TopographySurface
from Autodesk.Revit.UI.Selection import *
from System.Collections.Generic import *
from System.Collections import *
from System import *
from math import *

class ThisApplication (ApplicationEntryPoint):
	### CONFIGURATION #################################################################
	#
	# Metric or imperial?
	# True  -> Export unit meters
	# False -> Export unit inches
	metric = True
	#
	#
	# Output file normal Geometry
	fname = os.environ['USERPROFILE'] + r'\Desktop\revout.rad'
	#
	# Output file Topography Triangle Geometry
	topo_fname = os.environ['USERPROFILE'] + r'\Desktop\revout_topo.rad'
	#
	# Output file Topography Mesh Geometry
	# OBJ files with same base and numbered
	topo_obj_fname = os.environ['USERPROFILE'] + r'\Desktop\revout_topomesh.rad'
	#
	#
	### END CONFIGURATION #############################################################

	#region Revit Macros generated code
	def FinishInitialization(self):
		ApplicationEntryPoint.FinishInitialization(self)
		self.InternalStartup()
	def OnShutdown(self):
		self.InternalShutdown()
		ApplicationEntryPoint.OnShutdown(self)
	def InternalStartup(self):
		self.Startup()
	def InternalShutdown(self):
		self.Shutdown()
	#endregion
	def Startup(self): pass
	def Shutdown(self): pass

	def ExportToRadiance(self):
		if (self.ActiveUIDocument != None):
			self.__ExportToRadianceImplementation(self.Application, self.ActiveUIDocument.Document)
		else:
			self.__ExportToRadianceImplementation(self.Application, None)
	def __ExportToRadianceImplementation(self, app, doc):
		context  = _MyExportContext(doc, self.metric, self.fname)
		exporter = CustomExporter(doc, context)
		curview = doc.ActiveView
		exporter.Export(curview)

	def ExportTopoToRadiance(self):
		if (self.ActiveUIDocument != None):
			self.__ExportTopoToRadianceImplementation(self.Application, self.ActiveUIDocument.Document)
		else:
			self.__ExportTopoToRadianceImplementation(self.Application, None)
	def __ExportTopoToRadianceImplementation(self, app, doc):
		context  = _MyTopoExportContext(doc, self.metric, self.topo_fname)
		exporter = CustomExporter(doc, context)
		curview = doc.ActiveView
		exporter.Export(curview)

	def ExportTopoToObj(self):
		if (self.ActiveUIDocument != None):
			self.__ExportTopoToObjImplementation(self.Application, self.ActiveUIDocument.Document)
		else:
			self.__ExportTopoToObjImplementation(self.Application, None)
	def __ExportTopoToObjImplementation(self, app, doc):
		context  = _MyTopoObjExportContext(doc, self.metric, self.topo_obj_fname)
		exporter = CustomExporter(doc, context)
		curview = doc.ActiveView
		exporter.Export(curview)

	# Transaction mode
	def GetTransactionMode(self):
		return Attributes.TransactionMode.Manual
	# Addin Id
	def GetAddInId(self):
		return '98E9898F-B62F-415A-9AA3-1F705B44EAC0'


class _RadianceBaseExportContext(IExportContext):
	RND = 10

	def __init__(self, doc, metric, fname):
		self.fname = fname
		self.doc = doc
		self.cur_n = 0
		self.cancelled = False
		if metric:
			self.scale = 0.3048
			self.xforms = [Transform.Identity.ScaleBasis(self.scale)]
		else:
			self.scale = 1.0
			self.xforms = [Transform.Identity]

	def _norm_name(self, s):
		return s.replace(' ', '_')

	def swstr(self, lines):
		sl = ['# Software: %s %s' % (__software__, __version__),
			'# Created: %s' % time.asctime()]
		sl.extend(lines)
		sl.append('# Original File: %s\n\n' % self.doc.PathName)
		return '\n'.join(sl)

	def _show_exc(self):
		# Otherwise we'll only see the traceback to where C# called Python
		ei = sys.exc_info()
		sl = [str(ei[1])]
		tb = ei[2]
		while tb:
			s = '%s Line: %d' % (tb.tb_lasti, tb.tb_lineno)
			sl.append(s)
			tb = tb.tb_next
		ss = '\n'.join(sl)
		TaskDialog.Show('Traceback', ss)


	def _make_polygon(self, name, pts):
		sl = ['{name} polygon {name}_{n:06d}'.format(name=name, n=self.cur_n)]
		self.cur_n += 1
		sl.append('0 0 {num}'.format(num=len(pts)*3))
		for pt in pts:
			# The IronPython '{:g}'.format() is broken (uses locale),
			# so we need to use legacy string interpolation for floats
			sl.append('\t%.8g %.8g %.8g' % (round(pt[0],self.RND),
				round(pt[1],self.RND), round(pt[2],self.RND)))
		sl.append('')
		return '\n'.join(sl)

	def _make_ring(self, name, cen, norm, r0, r1):
		sl = ['{name} ring {name}_{n:06d}\n0 0 8'.format(
			name=name, n=self.cur_n)]
		self.cur_n += 1
		sl.append('\t%.8g %.8g %.8g' % (round(cen[0],self.RND),
			round(cen[1],self.RND), round(cen[2],self.RND)))
		sl.append('\t%.8g %.8g %.8g' % (round(norm[0],self.RND),
			round(norm[1],self.RND), round(norm[2],self.RND)))
		sl.append('\t%.8g %.8g\n' % (round(r0,self.RND), round(r1,self.RND)))
		return '\n'.join(sl)

	def _make_cylinder(self, name, start, end, rad):
		sl = ['{name} cylinder {name}_{n:06d}\n0 0 7'.format(
			name=name, n=self.cur_n)]
		self.cur_n += 1
		sl.append('\t%.8g %.8g %.8g' % (round(start[0],self.RND),
			round(start[1],self.RND), round(start[2],self.RND)))
		sl.append('\t%.8g %.8g %.8g' % (round(end[0],self.RND),
			round(end[1],self.RND), round(end[2],self.RND)))
		sl.append('\t%.8g\n' % round(rad,self.RND))
		return '\n'.join(sl)

	def _make_cone(self, name, start, end, r0, r1):
		sl = ['{name} cone {name}_{n:06d}\n0 0 8'.format(
			name=name, n=self.cur_n)]
		self.cur_n += 1
		sl.append('\t%.8g %.8g %.8g' % (round(start[0],self.RND),
			round(start[1],self.RND), round(start[2],self.RND)))
		sl.append('\t%.8g %.8g %.8g' % (round(end[0],self.RND),
			round(end[1],self.RND), round(end[2],self.RND)))
		sl.append('\t%.8g %.8g\n' % (round(r0,self.RND), round(r1,self.RND)))
		return '\n'.join(sl)

	def _make_objmesh(self, name, fname):
		s = '{name} mesh {name}_{n:06d}\n1 {fname}\n0\n0\n\n'.format(
				name=name, n=self.cur_n, fname=fname)
		self.cur_n += 1
		return s


class _MyTopoExportContext(_RadianceBaseExportContext):
	def __init__(self, doc, metric, fname):
		self.cur_mat = ''
		_RadianceBaseExportContext.__init__(self, doc, metric, fname)

	def Start(self):
		self.f = open(self.fname, 'w')
		self.f.write(self.swstr(['# Radiance Scene File',
			'# Revit Topography Mesh Export']))
		self.f.write('# Original File: %s\n\n' % self.doc.PathName)
		return True
	def Finish(self):
		self.f.write('\n# End of Revit Topography Mesh Export\n')
		self.f.close()
	def IsCanceled(self):
		if self.cancelled:
			self.f.write('\n# Export cancelled\n\n')
		return self.cancelled
	def OnDaylightPortal(self, node): pass
	def OnElementBegin(self, id):
		el = self.doc.GetElement(id)
		if isinstance(el, TopographySurface):
			matpar = el.Parameter['Material']
			if matpar and matpar.HasValue:
				matid = matpar.AsElementId()
				matel = self.doc.GetElement(matid)
				if matel:
					self.cur_mat = matel.Name
				elif el.Category and el.Category.Material: # assume by category
					self.cur_mat = el.Category.Material.Name
			return RenderNodeAction.Proceed
		return RenderNodeAction.Skip
	def OnElementEnd(self, node):
		self.cur_mat = ''
	def OnFaceBegin(self, node): return RenderNodeAction.Skip
	def OnFaceEnd(self, node): pass
	def OnInstanceBegin(self, node): return RenderNodeAction.Skip
	def OnInstanceEnd(self, node): pass
	def OnLight(self, node): pass
	def OnLinkBegin(self, node): return RenderNodeAction.Skip
	def OnLinkEnd(self, node): pass
	def OnMaterial(self, node): pass
	def OnPolymesh(self, node):
		xf = self.xforms[-1]
		if self.cur_mat:
			name = self._norm_name('mesh_topography+++' + self.cur_mat)
		else:
			name = 'mesh_topography'
		fs = node.GetFacets()
		points = node.GetPoints()
		self.f.write('\n# Polymesh Topography with %d facets\n'
				% node.NumberOfFacets)
		self.f.write('\n# points: %d - normals: %d - UVs: %d\n'
				% (node.NumberOfPoints, node.NumberOfNormals, node.NumberOfUVs))
		for facet in fs:
			pl = (xf.OfPoint(points[facet.V1]),
				xf.OfPoint(points[facet.V2]),
				xf.OfPoint(points[facet.V3]))
			self.f.write(self._make_polygon(name, pl))
	def OnRPC(self, node): pass
	def OnViewBegin(self, node): return RenderNodeAction.Proceed
	def OnViewEnd(self, node): pass

class _MyTopoObjExportContext(_RadianceBaseExportContext):
	def __init__(self, doc, metric, fname):
		self.cur_mat = ''
		rdot = fname.rfind('.')
		if rdot:
			self.basefn = fname[:rdot]
		else:
			self.basefn = fname
		_RadianceBaseExportContext.__init__(self, doc, metric, fname)
	def Start(self):
		self.f = open(self.fname, 'w')
		self.f.write(self.swstr(['# Radiance Mesh Files',
			'# Revit Topography Mesh Export']))
		return True
	def Finish(self):
		self.f.write('\n# End of Revit Topography Mesh Export\n')
		self.f.close()
	def IsCanceled(self):
		if self.cancelled:
			self.f.write('\n# Export cancelled\n\n')
		return self.cancelled
	def OnDaylightPortal(self, node): pass
	def OnElementBegin(self, id):
		el = self.doc.GetElement(id)
		if isinstance(el, TopographySurface):
			matpar = el.Parameter['Material']
			if matpar and matpar.HasValue:
				matid = matpar.AsElementId()
				matel = self.doc.GetElement(matid)
				if matel:
					self.cur_mat = matel.Name
				elif el.Category and el.Category.Material: # assume by category
					self.cur_mat = el.Category.Material.Name
			return RenderNodeAction.Proceed
		return RenderNodeAction.Skip
	def OnElementEnd(self, node):
		self.cur_mat = ''
	def OnFaceBegin(self, node): return RenderNodeAction.Skip
	def OnFaceEnd(self, node): pass
	def OnInstanceBegin(self, node): return RenderNodeAction.Skip
	def OnInstanceEnd(self, node): pass
	def OnLight(self, node): pass
	def OnLinkBegin(self, node): return RenderNodeAction.Skip
	def OnLinkEnd(self, node): pass
	def OnMaterial(self, node): pass
	def OnPolymesh(self, node):
		xf = self.xforms[-1]
		fullbase = '{base}_{n:06d}'.format(base=self.basefn, n=self.cur_n)
		if self.cur_mat:
			name = self._norm_name('mesh_topography+++' + self.cur_mat)
		else:
			name = 'mesh_topography'
		self.f.write(self._make_objmesh(name, fullbase + '.rfm'))

		objf = open(fullbase + '.obj', 'w')
		self.f.write(self.swstr(['# Wavefront Object File',
			'# Revit Topography Mesh Export for Radiance']))
		for v in node.GetPoints():
			v = xf.OfPoint(v)
			objf.write('v %g %g %g\n' % (round(v.X,self.RND),
				round(v.Y,self.RND), round(v.Z,self.RND)))
		objf.write('\n')
		for vt in node.GetUVs():
			objf.write('vt %g %g\n' % (round(vt.U*self.scale,self.RND),
				round(vt.V*self.scale,self.RND)))
		objf.write('\n')
		for vn in node.GetNormals():
			vn = xf.OfVector(vn)
			vn = vn.Normalize()
			objf.write('vn %g %g %g\n' % (round(vn.X,self.RND),
				round(vn.Y,self.RND), round(vn.Z,self.RND)))
		objf.write('\n')
		for f in node.GetFacets():
			objf.write(
					'f %(x)d/%(x)d/%(x)d %(y)d/%(y)d/%(y)d %(z)d/%(z)d/%(z)d\n'
					% {'x': f.V1+1, 'y': f.V2+1, 'z': f.V3+1})
		objf.write('\n')
		objf.write('# End of Revit Topography Mesh Export\n')
	def OnRPC(self, node): pass
	def OnViewBegin(self, node): return RenderNodeAction.Proceed
	def OnViewEnd(self, node): pass




class _MyExportContext(_RadianceBaseExportContext):
	def __init__(self, doc, metric, fname):
		self.cylhalfs = []
		self.circhalfs = []
		self.revhalfs = []
		self.cur_level = 'No-Level'
		self.cur_family = ['']
		self.cur_type = ''
		_RadianceBaseExportContext.__init__(self, doc, metric, fname)
	def Start(self):
		self.f = open(self.fname, 'w')
		self.f.write(self.swstr(['# Radiance Scene File',
			'# Revit Geometry Export']))
		return True
	def Finish(self):
		self.f.write('\n# End of Revit Geometry Export\n')
		self.f.close()
	def IsCanceled(self):
		if self.cancelled:
			self.f.write('\n# Export cancelled\n\n')
		return self.cancelled
	def OnDaylightPortal(self, node): pass

	def _get_element_level(self, el):
		level = self.doc.GetElement(el.LevelId)
		if level:
			return level.Name
		if hasattr(el, 'Host') and isinstance(el.Host, Level):
			return el.Host.Name
		# let's just hope this key is not locale dependent.
		ref_lev = el.Parameter['Reference Level']
		if ref_lev and ref_lev.HasValue:
			return self.doc.GetElement(ref_lev.AsElementId()).Name

	def _setup_halfcache(self):
		self.cylhalfs.append({})
		self.circhalfs.append({})
		self.revhalfs.append({})

	def _process_halfcache(self):
		# write the partial cylinders that were not matched
		cylhl = self.cylhalfs.pop()
		for hck, (ch, matid) in cylhl.items():
			self._write_cylinder_segments(ch, matid)
		# write the partial circles that were not matched
		circhl = self.circhalfs.pop()
		for hck, (ch, matid) in circhl.items():
			self._write_mesh_tris(ch.Triangulate(), matid)
		# write the partial cones that were not matched
		revhl = self.revhalfs.pop()
		for hck, (ch, matid) in revhl.items():
			self._write_mesh_tris(ch.Triangulate(), matid)

	def OnElementBegin(self, id):
		try:
			self._setup_halfcache()
			el = self.doc.GetElement(id)
			level = self._get_element_level(el)
			if level:
				self.cur_level = level
			else:
				self.cur_level = 'Level-XXX'

			cur_family = ''
			if hasattr(el, 'Symbol'):
				sym = el.Symbol
				if sym and hasattr(sym, 'Family') and sym.Family:
					cur_family = sym.Family.Name
					self.cur_family.append(cur_family)
			self.cur_type = el.Name
			self.f.write('\n# Element: %s // %s // %s\n'
					% (self.cur_level, self.cur_family, el.Name))
			return RenderNodeAction.Proceed
		except:
			self._show_exc()

	def OnElementEnd(self, id):
		self._process_halfcache()
		self.cur_level = 'No-Level'
		self.f.write('# On Element END: "%s"\n' % self.cur_family)
		if len(self.cur_family) > 1:
			self.cur_family.pop()

	def _make_name(self, matid):
		mat = self.doc.GetElement(matid)
		if self.cur_family[-1]:
			l = [self.cur_level, self.cur_family[-1]]
		else:
			l = [self.cur_level, self.cur_type]
		if mat:
			l.append(mat.Name)
		s = '+++'.join(l)
		return self._norm_name(s)

	def _write_mesh_tris(self, mesh, matid):
		mesh = mesh.Transformed[self.xforms[-1]]
		self.f.write('###### transforming in _write_mesh_tris()\n')
		name = self._make_name(matid)
		sl = []
		for i in range(mesh.NumTriangles):
			t = mesh.Triangle[i]
			s = self._make_polygon(name,
					(t.Vertex[0], t.Vertex[1], t.Vertex[2]))
			sl.append(s)
		self.f.write('\n### Mesh Face Triangles\n')
		self.f.write(''.join(sl))
		self.f.write('\n')

	def _almost_equal(self, x, y):
		return round(x, self.RND) == round(y, self.RND)

	def _loop_to_segments(self, loop):
		segs = []
		for edge in loop:
			c = edge.AsCurve().CreateTransformed(self.xforms[-1])
			cpts = c.Tessellate()
			segs.append(list(cpts))
		self.f.write('###### transforming in _loop_to_segments()\n')
		for i in range(len(segs)-1):
			# some edges may be backwards
			pl0 = segs[i]
			pl1 = segs[i+1]
			if pl0[-1].IsAlmostEqualTo( pl1[0]):
				continue
			if pl0[-1].IsAlmostEqualTo(pl1[-1]):
				pl1.reverse()
				continue
			if pl0[0].IsAlmostEqualTo(pl1[0]):
				pl0.reverse()
				continue
			if pl0[0].IsAlmostEqualTo(pl1[-1]):
				pl0.reverse()
				pl1.reverse()
				continue
		return segs

	def _write_polygon(self, f, matid):
		name = self._make_name(matid)
		edges = f.EdgeLoops[0]
		loop = edges.ForwardIterator()
		segs = self._loop_to_segments(loop)
		pts = []
		for ptl in segs:
			pts.extend(ptl[:-1])
		if f.EdgeLoops.Size > 1: # holes
			last = pts[-1]
			for i in range(1, f.EdgeLoops.Size):
				loop = f.EdgeLoops[i]
				hpts = []
				hsegs = self._loop_to_segments(loop)
				for hptl in hsegs:
					hpts.extend(hptl[:-1])
				hpts.append(hpts[0])
				hpts.reverse()
				hpts.append(last)
				pts.extend(hpts)
		self.f.write(self._make_polygon(name, pts))

	def _write_circle(self, f, matid):
		name = self._make_name(matid)
		xf = self.xforms[-1]
		loop = f.EdgeLoops[0]
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			if isinstance(c, Arc):
				break
		curve = c.CreateTransformed(xf)
		self.f.write('###### transforming in write_circle()\n')
		irad = 0.0
		if f.EdgeLoops.Size == 2:
			curve2 = f.EdgeLoops[1][0].AsCurve().CreateTransformed(xf)
			irad = curve2.Radius
		s = self._make_ring(name, curve.Center, curve.Normal.Normalize(),
				irad, curve.Radius)
		self.f.write(s)

	def _is_circle(self, f):
		loops = f.EdgeLoops
		if loops.Size not in [1,2]:
			return False
		loop = loops[0]
		if loop.Size != 2:
			return False
		c0 = loop[0].AsCurve()
		c1 = loop[1].AsCurve()
		if not isinstance(c0, Arc) or not isinstance(c1, Arc):
			return False
		if not self._almost_equal(c0.Radius, c1.Radius):
			return False
		if loops.Size == 2: # inner loop for ring
			loop2 = loops[1]
			if loop2.Size != 2:
				return False
			c02 = loop2[0].AsCurve()
			c12 = loop2[1].AsCurve()
			if not isinstance(c02, Arc) or not isinstance(c12, Arc):
				return False
			if not self._almost_equal(c02.Radius, c12.Radius):
				return False
			if not c0.Center.IsAlmostEqualTo(c02.Center):
				return False
		return True

	def _get_hcirc_sig(self, cyl):
		# unique signature for the location and size of a circle
		# we use this to reunite half shells
		# XXX Should we transform them to world first ?!?
		# XXX Since this works on a per-family instance basis,
		# XXX that may not actually be necessary.
		loop = cyl.EdgeLoops[0]
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			if isinstance(c, Arc):
				break
		orig = c.Center
		norm = c.Center
		rad = c.Radius
		return (round(rad, self.RND),
			round(orig[0], self.RND),
			round(orig[1], self.RND),
			round(orig[2], self.RND),
			round(norm[0], self.RND),
			round(norm[1], self.RND),
			round(norm[2], self.RND))

	def _is_half_circle(self,f):
		loop = f.EdgeLoops[0]
		if loop.Size != 3:
			return False
		curves = []
		lines = []
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			if isinstance(c, Arc):
				curves.append(c)
			elif isinstance(c, Line):
				lines.append(c)
			else:
				return False
		if (not len(curves) == 1) or (not len(lines) == 2):
			return False
		c = curves[0]
		l = lines[0]
		if ((not c.Center.IsAlmostEqualTo(l.GetEndPoint(0)))
				and (not c.Center.IsAlmostEqualTo(l.GetEndPoint(1)))):
			return
		c = curves[0]
		if (round(c.GetEndParameter(1) - c.GetEndParameter(0), self.RND) * 2
				== round(c.Period, self.RND)):
			return True
		return False

	def _is_half_cylinder(self,f):
		loop = f.EdgeLoops[0]
		if loop.Size != 4:
			self.f.write('\n# XXX Cylinder with weird number of edges #%d\n'
					% loop.Size)
			return 'weird'
		curves = []
		lines = []
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			if isinstance(c, Arc):
				curves.append(c)
			elif isinstance(c, Line):
				lines.append(c)
			else:
				self.f.write(
						'\n# XXX Cylinder with curved edge type #%d: "%s"\n'
						% (i, str(c.GetType())))
				return 'weird'
		if not len(curves) == len(lines) == 2:
			return False
		c = curves[0]
		if (abs(round(c.GetEndParameter(1) - c.GetEndParameter(0), self.RND)*2)
				== round(c.Period, self.RND)):
			return True
		return False

	def _is_half_rev(self,f):
		loop = f.EdgeLoops[0]
		if loop.Size not in [3,4]:
			self.f.write('# Rev loops size: %d\n' % loop.Size)
			return False
		arcs = []
		lines = []
		arc_i = -1
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			fullangle = c.GetEndParameter(1) - c.GetEndParameter(0)
			if (isinstance(c, Arc)
					and round(fullangle,self.RND)* 2 == round(c.Period,self.RND)
					and c.Normal.CrossProduct(f.Axis).IsZeroLength()):
				arcs.append(c)
			else:
				lines.append(c)
		if (not len(arcs) in [1,2]) or (not len(lines) == 2):
			self.f.write('# Rev arcs: %d - lines: %d\n'
					% (len(arcs), len(lines)))
			return False
		return True

	def _get_hrev_sig(self, cyl):
		# unique signature for the location and size of a revolved cone
		# we use this to reunite half shells
		# XXX Should we transform them to world first ?!?
		# XXX Since this works on a per-family instance basis,
		# XXX that may not actually be necessary.
		loop = cyl.EdgeLoops[0]
		cl = []
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			fullangle = c.GetEndParameter(1) - c.GetEndParameter(0)
			if (isinstance(c, Arc)
					and round(fullangle,self.RND)*2 ==round(c.Period,self.RND)):
				cl.append(c)
				break
		c0 = cl[0].Center
		if len(cl) == 2:
			c1 = cl[1].Center
			r1 = cl[1].Radius
		else:
			c1 = (0.0,0.0,0.0)
			r1 = 0.0
		return (round(cl[0].Radius, self.RND),
			round(r1, self.RND),
			round(c0[0], self.RND),
			round(c0[1], self.RND),
			round(c0[2], self.RND),
			round(c1[0], self.RND),
			round(c1[1], self.RND),
			round(c1[2], self.RND))

	def _write_rev(self, face, matid):
		name = self._make_name(matid)
		loop = face.EdgeLoops[0]
		arcs = []
		lines = []
		arc_i = -1
		for i in range(loop.Size):
			c = loop[i].AsCurve()
			fullangle = c.GetEndParameter(1) - c.GetEndParameter(0)
			if (isinstance(c, Arc)
					and round(fullangle,self.RND)*2 ==round(c.Period,self.RND)):
				arcs.append(c)
				if arc_i == -1:
					arc_i = i
			else:
				lines.append(c)
		xf = self.xforms[-1]
		segs = self._loop_to_segments(loop)
		if arc_i == 0:
			seg0 = segs[1]
			seg1 = segs[-1]
		else:
			seg0 = segs[0]
			seg1 = segs[2]
		seg1.reverse()
		rl = [p0.DistanceTo(p1)*.5 for p0, p1 in zip(seg0, seg1)]
		cl = [p0.Add(p1).Divide(2.0) for p0, p1 in zip(seg0, seg1)]
		self.f.write('###### transforming center points in _write_rev()\n')
		for i in range(len(rl)-1):
			c0 = xf.OfPoint(cl[i])
			c1 = xf.OfPoint(cl[i+1])
			r0 = rl[i]
			r1 = rl[i+1]
			s = self._make_cone(name, start=cl[i], end=cl[i+1], r0=r0, r1=r1)
			self.f.write(s)

	def _write_cylinder(self, cyl, matid):
		name = self._make_name(matid)
		arcs = self._get_cyl_arcs(cyl)
		xf = self.xforms[-1]
		self.f.write('###### transforming in _write_cylinder()\n')
		s = self._make_cylinder(name=self._make_name(matid),
			start=xf.OfPoint(arcs[0].Center),
			end=xf.OfPoint(arcs[1].Center),
			rad=xf.OfVector(cyl.Radius[0]).GetLength())
		self.f.write(s)

	def _get_cyl_arcs(self, cyl):
		loop = cyl.EdgeLoops[0]
		c0 = loop[0].AsCurve()
		if isinstance(c0, Arc):
			c1 = loop[2].AsCurve()
		else:
			c0 = loop[1].AsCurve()
			assert(isinstance(c0, Arc))
			c1 = loop[3].AsCurve()
		return c0, c1

	def _get_cyl_segments(self, cyl):
		loop = cyl.EdgeLoops[0]
		segs = self._loop_to_segments(loop)
		c0 = loop[0].AsCurve()
		if isinstance(c0, Arc):
			s0 = segs[0]
			s1 = segs[2]
		else:
			s0 = segs[1]
			s1 = segs[3]
		return s0, s1

	def _write_cylinder_segments(self, cyl, matid):
		name = self._make_name(matid)
		segs = self._get_cyl_segments(cyl)
		pts0 = segs[0]
		pts1 = segs[1]
		sl = []
		llen = pts0.Count
		if pts0.Count != pts1.Count:
			llen = min((pts0.Count, pts1.Count))
			self.f.write('# XXX Segment count differs: %d %d\n'
					% (pts0.Count, pts1.Count))
		for i in range(llen-1):
			s = self._make_polygon(name,
				(pts0[i], pts0[i+1], pts1[llen-i-2], pts1[llen-i-1]))
			sl.append(s)
			self.cur_n += 1
		self.f.write('\n### Cylinder Segments\n')
		self.f.write(''.join(sl))
		self.f.write('\n')

	def _get_hcyl_sig(self, cyl):
		# unique signature for the location and size of a cylinder
		# we use this to reunite half shells
		# XXX Should we transform them to world first ?!?
		# XXX Since this works on a per-family instance basis,
		# XXX that may not actually be necessary.
		arcs = self._get_cyl_arcs(cyl)
		orig = arcs[0].Center
		end = arcs[1].Center
		rad = cyl.Radius[0].GetLength()
		return (round(rad, self.RND),
			round(orig[0], self.RND),
			round(orig[1], self.RND),
			round(orig[2], self.RND),
			round(end[0], self.RND),
			round(end[1], self.RND),
			round(end[2], self.RND))

	def _process_cylinder(self, f, matid):
		hcres = self._is_half_cylinder(f)
		if hcres is True:
			hcsig = self._get_hcyl_sig(f)
			if self.cylhalfs and self.cylhalfs[-1].get(hcsig):
				del self.cylhalfs[-1][hcsig]
				self._write_cylinder(f, matid)
			elif self.cylhalfs:
				self.cylhalfs[-1][hcsig] = (f, matid)
			else:
				self._write_cylinder_segments(f, matid)
		elif hcres == 'weird':
			# unexpected properties for a cylyndrical face
			# probably a curved cylinder ("flexible tube")
			self._write_mesh_tris(f.Triangulate(), matid)
		else:
			self._write_cylinder_segments(f, matid)

	def _process_planar(self, f, matid):
		if self._is_circle(f):
			self._write_circle(f, matid)
		elif self._is_half_circle(f):
			hcsig = self._get_hcirc_sig(f)
			oh = self.circhalfs[-1].get(hcsig)
			if oh:
				del self.circhalfs[-1][hcsig]
				self._write_circle(f, matid)
			else:
				self.circhalfs[-1][hcsig] = (f, matid)
		else:
			self._write_polygon(f, matid)

	def _process_rev(self, f, matid):
		hcres = self._is_half_rev(f)
		if hcres is True:
			hcsig = self._get_hrev_sig(f)
			oh = self.revhalfs[-1].get(hcsig)
			if oh:
				del self.revhalfs[-1][hcsig]
				self._write_rev(f, matid)
			else:
				self.revhalfs[-1][hcsig] = (f, matid)
		else:
			self._write_mesh_tris(f.Triangulate(), matid)

	def OnFaceBegin(self, node):
		try:
			f = node.GetFace()
			ref = f.Reference
			matid = f.MaterialElementId
			if isinstance(f, CylindricalFace):
				self._process_cylinder(f, matid)
			elif isinstance(f, PlanarFace):
				self._process_planar(f, matid)
			elif isinstance(f, RevolvedFace):
				self._process_rev(f, matid)
			elif isinstance(f, ConicalFace):
				self._process_rev(f, matid)
			else:
				self.f.write('### Face type: %s\n' % str(f.GetType()))
				self._write_mesh_tris(f.Triangulate(), matid)
			return RenderNodeAction.Skip
		except:
			self._show_exc()


	def OnFaceEnd(self, node): pass

	def OnInstanceBegin(self, node):
		sym = self.doc.GetElement(node.GetSymbolId())
		if sym and hasattr(sym, 'Family') and sym.Family:
			cur_family = sym.Family.Name
			self.f.write('###---### Instance of: "%s"\n' % cur_family)
		self.xforms.append(self.xforms[-1].Multiply(node.GetTransform()))
		self._setup_halfcache()
		return RenderNodeAction.Proceed
	def OnInstanceEnd(self, node):
		self._process_halfcache()
		sym = self.doc.GetElement(node.GetSymbolId())
		if sym and hasattr(sym, 'Family') and sym.Family:
			cur_family = sym.Family.Name
			self.f.write('###---### END Instance of: "%s"\n' % cur_family)
		self.xforms.pop()


	def OnLight(self, node): pass
	def OnLinkBegin(self, node):
		TaskDialog.Show('OnLinkBegin', str(node))
		return RenderNodeAction.Proceed
	def OnLinkEnd(self, node): pass
	def OnMaterial(self, node): pass
	def OnPolymesh(self, node): pass # seperate function above
	def OnRPC(self, node): pass
	def OnViewBegin(self, node): return RenderNodeAction.Proceed
	def OnViewEnd(self, node): pass


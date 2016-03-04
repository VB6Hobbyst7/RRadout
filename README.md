# RRadout
RRadout is a Python module for Autodesk Revit, exporting geometry data to the 
<a href="http://www.radiance-online.org/">Radiance</a> lighting simulation package.

<b>Warning:</b> This is beta software! (Release 0.2 has advanced from "Proof of Concept" to "Prototype" status.)

Basic features of RRadout are:

<ul TYPE=SQUARE>
	<li> Export all normal geometry visible in the active view (except RPCs)
	</li>
	<li> Circular PlanarFaces as rings where possible (including coaxial hole).
	</li>
	<li> Other PlanarFaces as polygons (including holes)
	</li>
	<li> CylindricalFaces as cylinders where possible (otherwise rectangles).
	</li>
	<li> ConicalFaces as cones where possible (otherwise triangles).
	</li>
	<li> RevolvedFaces as (sequence of) cones where possible (otherwise
		triangles).
	</li>
	<li> All other surfaces tessellated into triangles.
	</li>
	<li> Names based on level, leaf node family name (type name for
		non-family types), and material
	</li>
	<li> Export in metric or imperial units (configured in code)
	</li>
	<li> Separate function for exporting topography meshes as triangles.
	</li>
	<li> Separate function for exporting topography meshes as Wavefront
		Object files for use with obj2mesh (with normals for smoothing and
		UV-coordinates scaled to model units).
	</li>
</ul>

More details on the <a href="http://www.schorsch.com/en/download/rradout/">RRadout Homepage</a>

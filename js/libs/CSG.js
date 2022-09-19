// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
// 
// Example usage:
// 
//     var cube = CSG.cube();
//     var sphere = CSG.sphere({ radius: 1.3 });
//     var polygons = cube.subtract(sphere).toPolygons();
// 
// ## Implementation Details
// 
// All CSG operations are implemented in terms of two functions, `clipTo()` and
// `invert()`, which remove parts of a BSP tree inside another BSP tree and swap
// solid and empty space, respectively. To find the union of `a` and `b`, we
// want to remove everything in `a` inside `b` and everything in `b` inside `a`,
// then combine polygons from `a` and `b` into one solid:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     a.build(b.allPolygons());
// 
// The only tricky part is handling overlapping coplanar polygons in both trees.
// The code above keeps both copies, but we need to keep them in one tree and
// remove them in the other tree. To remove them from `b` we can clip the
// inverse of `b` against `a`. The code for union now looks like this:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     b.invert();
//     b.clipTo(a);
//     b.invert();
//     a.build(b.allPolygons());
// 
// Subtraction and intersection naturally follow from set operations. If
// union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
// `A & B = ~(~A | ~B)` where `~` is the complement operator.
// 
// ## License
// 
// Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the MIT license.

// # class CSG

// Holds a binary space partition tree representing a 3D solid. Two solids can
// be combined using the `union()`, `subtract()`, and `intersect()` methods.

CSG = function () {
  this.polygons = [];
};

class Vector {
  constructor(x=0, y=0, z=0) {
      this.x=x;
      this.y=y;
      this.z=z;
  }
  copy(v){
      this.x=v.x;
      this.y=v.y;
      this.z=v.z;
      return this
  }
  clone() {
      return new Vector(this.x,this.y,this.z)
  }
  negate() {
      this.x*=-1;
      this.y*=-1;
      this.z*=-1;
      return this
  }
  add(a) {
      this.x+=a.x
      this.y+=a.y
      this.z+=a.z
      return this;
  }
  sub(a) {
      this.x-=a.x
      this.y-=a.y
      this.z-=a.z
      return this
  }
  times(a) {
      this.x*=a
      this.y*=a
      this.z*=a
      return this
  }
  dividedBy(a) {
      this.x/=a
      this.y/=a
      this.z/=a
      return this
  }
  lerp(a, t) {
      return this.add(tv0.copy(a).sub(this).times(t))
  }
  unit() {
      return this.dividedBy(this.length())
  }
  length(){
      return Math.sqrt((this.x**2)+(this.y**2)+(this.z**2))
  }
  normalize(){
      return this.unit()
  }
  cross(b) {
      let a = this;
  const ax = a.x, ay = a.y, az = a.z;
  const bx = b.x, by = b.y, bz = b.z;

  this.x = ay * bz - az * by;
  this.y = az * bx - ax * bz;
  this.z = ax * by - ay * bx;

  return this;
  }
  dot(b){
      return (this.x*b.x)+(this.y*b.y)+(this.z*b.z)
  }
}

// Construct a CSG solid from a list of `CSG.Polygon` instances.
CSG.fromPolygons = function (polygons) {
  var csg = new CSG();
  csg.polygons = polygons;
  return csg;
};

CSG.prototype = {
  clone: function () {
    var csg = new CSG();
    csg.polygons = this.polygons.map(function (p) { return p.clone(); });
    return csg;
  },

  toPolygons: function () {
    return this.polygons;
  },

  // Return a new CSG solid representing space in either this solid or in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.union(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |       +----+
  //     +----+--+    |       +----+       |
  //          |   B   |            |       |
  //          |       |            |       |
  //          +-------+            +-------+
  // 
  union: function (csg) {
    var a = new CSG.Node(this.clone().polygons);
    var b = new CSG.Node(csg.clone().polygons);
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.build(b.allPolygons());
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid representing space in this solid but not in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.subtract(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |    +--+
  //     +----+--+    |       +----+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  subtract: function (csg) {
    var a = new CSG.Node(this.clone().polygons);
    var b = new CSG.Node(csg.clone().polygons);
    a.invert();
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.build(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid representing space both this solid and in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.intersect(B)
  // 
  //     +-------+
  //     |       |
  //     |   A   |
  //     |    +--+----+   =   +--+
  //     +----+--+    |       +--+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  intersect: function (csg) {
    var a = new CSG.Node(this.clone().polygons);
    var b = new CSG.Node(csg.clone().polygons);
    a.invert();
    b.clipTo(a);
    b.invert();
    a.clipTo(b);
    b.clipTo(a);
    a.build(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid with solid and empty space switched. This solid is
  // not modified.
  inverse: function () {
    var csg = this.clone();
    csg.polygons.map(function (p) { p.flip(); });
    return csg;
  }
};

// Construct an axis-aligned solid cuboid. Optional parameters are `center` and
// `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
// specified using a single number or a list of three numbers, one for each axis.
// 
// Example code:
// 
//     var cube = CSG.cube({
//       center: [0, 0, 0],
//       radius: 1
//     });
CSG.cube = function (options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = !options.radius ? [1, 1, 1] : options.radius.length ?
    options.radius : [options.radius, options.radius, options.radius];
  return CSG.fromPolygons([
    [[0, 4, 6, 2], [-1, 0, 0]],
    [[1, 3, 7, 5], [+1, 0, 0]],
    [[0, 1, 5, 4], [0, -1, 0]],
    [[2, 6, 7, 3], [0, +1, 0]],
    [[0, 2, 3, 1], [0, 0, -1]],
    [[4, 5, 7, 6], [0, 0, +1]]
  ].map(function (info) {
    return new CSG.Polygon(info[0].map(function (i) {
      var pos = new CSG.Vector(
        c.x + r[0] * (2 * !!(i & 1) - 1),
        c.y + r[1] * (2 * !!(i & 2) - 1),
        c.z + r[2] * (2 * !!(i & 4) - 1)
      );
      return new CSG.Vertex(pos, new CSG.Vector(info[1]));
    }));
  }));
};

// Construct a solid sphere. Optional parameters are `center`, `radius`,
// `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
// The `slices` and `stacks` parameters control the tessellation along the
// longitude and latitude directions.
// 
// Example usage:
// 
//     var sphere = CSG.sphere({
//       center: [0, 0, 0],
//       radius: 1,
//       slices: 16,
//       stacks: 8
//     });
CSG.sphere = function (options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var stacks = options.stacks || 8;
  var polygons = [], vertices;
  function vertex(theta, phi) {
    theta *= Math.PI * 2;
    phi *= Math.PI;
    var dir = new CSG.Vector(
      Math.cos(theta) * Math.sin(phi),
      Math.cos(phi),
      Math.sin(theta) * Math.sin(phi)
    );
    vertices.push(new CSG.Vertex(c.plus(dir.times(r)), dir));
  }
  for (var i = 0; i < slices; i++) {
    for (var j = 0; j < stacks; j++) {
      vertices = [];
      vertex(i / slices, j / stacks);
      if (j > 0) vertex((i + 1) / slices, j / stacks);
      if (j < stacks - 1) vertex((i + 1) / slices, (j + 1) / stacks);
      vertex(i / slices, (j + 1) / stacks);
      polygons.push(new CSG.Polygon(vertices));
    }
  }
  return CSG.fromPolygons(polygons);
};

// Construct a solid cylinder. Optional parameters are `start`, `end`,
// `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
// `16`. The `slices` parameter controls the tessellation.
// 
// Example usage:
// 
//     var cylinder = CSG.cylinder({
//       start: [0, -1, 0],
//       end: [0, 1, 0],
//       radius: 1,
//       slices: 16
//     });
CSG.cylinder = function (options) {
  options = options || {};
  var s = new CSG.Vector(options.start || [0, -1, 0]);
  var e = new CSG.Vector(options.end || [0, 1, 0]);
  var ray = e.minus(s);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var axisZ = ray.unit(), isY = (Math.abs(axisZ.y) > 0.5);
  var axisX = new CSG.Vector(isY, !isY, 0).cross(axisZ).unit();
  var axisY = axisX.cross(axisZ).unit();
  var start = new CSG.Vertex(s, axisZ.negated());
  var end = new CSG.Vertex(e, axisZ.unit());
  var polygons = [];
  function point(stack, slice, normalBlend) {
    var angle = slice * Math.PI * 2;
    var out = axisX.times(Math.cos(angle)).plus(axisY.times(Math.sin(angle)));
    var pos = s.plus(ray.times(stack)).plus(out.times(r));
    var normal = out.times(1 - Math.abs(normalBlend)).plus(axisZ.times(normalBlend));
    return new CSG.Vertex(pos, normal);
  }
  for (var i = 0; i < slices; i++) {
    var t0 = i / slices, t1 = (i + 1) / slices;
    polygons.push(new CSG.Polygon([start, point(0, t0, -1), point(0, t1, -1)]));
    polygons.push(new CSG.Polygon([point(0, t1, 0), point(0, t0, 0), point(1, t0, 0), point(1, t1, 0)]));
    polygons.push(new CSG.Polygon([end, point(1, t1, 1), point(1, t0, 1)]));
  }
  return CSG.fromPolygons(polygons);
};

// # class Vector

// Represents a 3D vector.
// 
// Example usage:
// 
//     new CSG.Vector(1, 2, 3);
//     new CSG.Vector([1, 2, 3]);
//     new CSG.Vector({ x: 1, y: 2, z: 3 });

CSG.Vector = function (x, y, z) {
  if (arguments.length == 3) {
    this.x = x;
    this.y = y;
    this.z = z;
  } else if ('x' in x) {
    this.x = x.x;
    this.y = x.y;
    this.z = x.z;
  } else {
    this.x = x[0];
    this.y = x[1];
    this.z = x[2];
  }
};

CSG.Vector.prototype = {
  clone: function () {
    return new CSG.Vector(this.x, this.y, this.z);
  },

  negated: function () {
    return new CSG.Vector(-this.x, -this.y, -this.z);
  },

  plus: function (a) {
    return new CSG.Vector(this.x + a.x, this.y + a.y, this.z + a.z);
  },

  minus: function (a) {
    return new CSG.Vector(this.x - a.x, this.y - a.y, this.z - a.z);
  },

  times: function (a) {
    return new CSG.Vector(this.x * a, this.y * a, this.z * a);
  },

  dividedBy: function (a) {
    return new CSG.Vector(this.x / a, this.y / a, this.z / a);
  },

  dot: function (a) {
    return this.x * a.x + this.y * a.y + this.z * a.z;
  },

  lerp: function (a, t) {
    return this.plus(a.minus(this).times(t));
  },

  length: function () {
    return Math.sqrt(this.dot(this));
  },

  unit: function () {
    return this.dividedBy(this.length());
  },

  cross: function (a) {
    return new CSG.Vector(
      this.y * a.z - this.z * a.y,
      this.z * a.x - this.x * a.z,
      this.x * a.y - this.y * a.x
    );
  }
};

// # class Vertex

// Represents a vertex of a polygon. Use your own vertex class instead of this
// one to provide additional features like texture coordinates and vertex
// colors. Custom vertex classes need to provide a `pos` property and `clone()`,
// `flip()`, and `interpolate()` methods that behave analogous to the ones
// defined by `CSG.Vertex`. This class provides `normal` so convenience
// functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
// is not used anywhere else.

CSG.Vertex = function (pos, normal) {
  this.pos = new CSG.Vector(pos);
  this.normal = new CSG.Vector(normal);
};

CSG.Vertex.prototype = {
  clone: function () {
    return new CSG.Vertex(this.pos.clone(), this.normal.clone());
  },

  // Invert all orientation-specific data (e.g. vertex normal). Called when the
  // orientation of a polygon is flipped.
  flip: function () {
    this.normal = this.normal.negated();
  },

  // Create a new vertex between this vertex and `other` by linearly
  // interpolating all properties using a parameter of `t`. Subclasses should
  // override this to interpolate additional properties.
  interpolate: function (other, t) {
    return new CSG.Vertex(
      this.pos.lerp(other.pos, t),
      this.normal.lerp(other.normal, t)
    );
  }
};

// # class Plane

// Represents a plane in 3D space.

CSG.Plane = function (normal, w) {
  this.normal = normal;
  this.w = w;
};

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
CSG.Plane.EPSILON = 1e-5;

CSG.Plane.fromPoints = function (a, b, c) {
  var n = b.minus(a).cross(c.minus(a)).unit();
  return new CSG.Plane(n, n.dot(a));
};

CSG.Plane.prototype = {
  clone: function () {
    return new CSG.Plane(this.normal.clone(), this.w);
  },

  flip: function () {
    this.normal = this.normal.negated();
    this.w = -this.w;
  },

  // Split `polygon` by this plane if needed, then put the polygon or polygon
  // fragments in the appropriate lists. Coplanar polygons go into either
  // `coplanarFront` or `coplanarBack` depending on their orientation with
  // respect to this plane. Polygons in front or in back of this plane go into
  // either `front` or `back`.
  splitPolygon: function (polygon, coplanarFront, coplanarBack, front, back) {
    var COPLANAR = 0;
    var FRONT = 1;
    var BACK = 2;
    var SPANNING = 3;

    // Classify each point as well as the entire polygon into one of the above
    // four classes.
    var polygonType = 0;
    var types = [];
    for (var i = 0; i < polygon.vertices.length; i++) {
      var t = this.normal.dot(polygon.vertices[i].pos) - this.w;
      var type = (t < -CSG.Plane.EPSILON) ? BACK : (t > CSG.Plane.EPSILON) ? FRONT : COPLANAR;
      polygonType |= type;
      types.push(type);
    }

    // Put the polygon in the correct list, splitting it when necessary.
    switch (polygonType) {
      case COPLANAR:
        (this.normal.dot(polygon.plane.normal) > 0 ? coplanarFront : coplanarBack).push(polygon);
        break;
      case FRONT:
        front.push(polygon);
        break;
      case BACK:
        back.push(polygon);
        break;
      case SPANNING:
        var f = [], b = [];
        for (var i = 0; i < polygon.vertices.length; i++) {
          var j = (i + 1) % polygon.vertices.length;
          var ti = types[i], tj = types[j];
          var vi = polygon.vertices[i], vj = polygon.vertices[j];
          if (ti != BACK) f.push(vi);
          if (ti != FRONT) b.push(ti != BACK ? vi.clone() : vi);
          if ((ti | tj) == SPANNING) {
            var t = (this.w - this.normal.dot(vi.pos)) / this.normal.dot(vj.pos.minus(vi.pos));
            var v = vi.interpolate(vj, t);
            f.push(v);
            b.push(v.clone());
          }
        }
        if (f.length >= 3) front.push(new CSG.Polygon(f, polygon.shared));
        if (b.length >= 3) back.push(new CSG.Polygon(b, polygon.shared));
        break;
    }
  }
};

// # class Polygon

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).

CSG.Polygon = function (vertices) {
  this.vertices = vertices;
  this.plane = CSG.Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos);
};

CSG.Polygon.prototype = {
  clone: function () {
    var vertices = this.vertices.map(function (v) { return v.clone(); });
    return new CSG.Polygon(vertices);
  },

  flip: function () {
    this.vertices.reverse().map(function (v) { v.flip(); });
    this.plane.flip();
  }
};

let nbuf3=(ct)=>{
  return{
      top:0,
      array:new Float32Array(ct),
      write:function(v){(this.array[this.top++]=v.x);(this.array[this.top++]=v.y);(this.array[this.top++]=v.z);}
  }
}

CSG.toGeometry = function (csg, buffered = true) {
  let ps = csg.polygons;
  let geom;
  let g2;
  if (!buffered) //Old geometry path...
  {
    geom = new Geometry();
    let vs = geom.vertices;
    let fvuv = geom.faceVertexUvs[0]
    for (let i = 0; i < ps.length; i++) {
      let p = ps[i]
      let pvs = p.vertices;
      let v0 = vs.length;
      let pvlen = pvs.length

      for (let j = 0; j < pvlen; j++)
        vs.push(new THREE.Vector3().copy(pvs[j].pos))

      for (let j = 3; j <= pvlen; j++) {
        let fc = new THREE.Face3();
        let fuv = []
        fvuv.push(fuv)
        let fnml = fc.vertexNormals;
        fc.a = v0;
        fc.b = v0 + j - 2;
        fc.c = v0 + j - 1;

        fnml.push(new THREE.Vector3().copy(pvs[0].normal))
        fnml.push(new THREE.Vector3().copy(pvs[j - 2].normal))
        fnml.push(new THREE.Vector3().copy(pvs[j - 1].normal))
        fuv.push(new THREE.Vector3().copy(pvs[0].uv))
        fuv.push(new THREE.Vector3().copy(pvs[j - 2].uv))
        fuv.push(new THREE.Vector3().copy(pvs[j - 1].uv))

        fc.normal = new THREE.Vector3().copy(p.plane.normal)
        geom.faces.push(fc)
      }
    }
    geom = new THREE.BufferGeometry().fromGeometry(geom)
    geom.verticesNeedUpdate = geom.elementsNeedUpdate = geom.normalsNeedUpdate = true;
  } else { //BufferGeometry path
    let triCount = 0;
    ps.forEach(p => triCount += (p.vertices.length - 2))
    geom = new THREE.BufferGeometry()

    let vertices = nbuf3(triCount * 3 * 3)
    let normals = nbuf3(triCount * 3 * 3)
    let uvs; // = nbuf2(triCount * 2 * 3)
    let colors;
    let grps = []
    ps.forEach(p => {
      let pvs = p.vertices
      let pvlen = pvs.length
      if (p.shared !== undefined) {
        if (!grps[p.shared]) grps[p.shared] = []
      }
      if (pvlen) {
        if (pvs[0].color !== undefined) {
          if (!colors) colors = nbuf3(triCount * 3 * 3);
        }
        if (pvs[0].uv !== undefined) {
          if (!uvs) uvs = nbuf2(triCount * 2 * 3)
        }
      }
      for (let j = 3; j <= pvlen; j++) {
        (p.shared !== undefined) && (grps[p.shared].push(vertices.top / 3, (vertices.top / 3) + 1, (vertices.top / 3) + 2));
        vertices.write(pvs[0].pos)
        vertices.write(pvs[j - 2].pos)
        vertices.write(pvs[j - 1].pos)
        normals.write(pvs[0].normal)
        normals.write(pvs[j - 2].normal)
        normals.write(pvs[j - 1].normal);
        uvs && (pvs[0].uv) && (uvs.write(pvs[0].uv) || uvs.write(pvs[j - 2].uv) || uvs.write(pvs[j - 1].uv));
        colors && (colors.write(pvs[0].color) || colors.write(pvs[j - 2].color) || colors.write(pvs[j - 1].color))
      }
    }
    )
    geom.setAttribute('position', new THREE.BufferAttribute(vertices.array, 3));
    geom.setAttribute('normal', new THREE.BufferAttribute(normals.array, 3));
    uvs && geom.setAttribute('uv', new THREE.BufferAttribute(uvs.array, 2));
    colors && geom.setAttribute('color', new THREE.BufferAttribute(colors.array, 3));
    if (grps.length) {
      let index = []
      let gbase = 0;
      for (let gi = 0; gi < grps.length; gi++) {
        geom.addGroup(gbase, grps[gi].length, gi)
        gbase += grps[gi].length
        index = index.concat(grps[gi]);
      }
      geom.setIndex(index)
    }
    g2 = geom;
  }
  return geom;
}

CSG.fromGeometry = function(geom) {
  let polys = []
  if (geom.isGeometry) {
      let fs = geom.faces;
      let vs = geom.vertices;
      let fm = ['a', 'b', 'c']
      for (let i = 0; i < fs.length; i++) {
          let f = fs[i];
          let vertices = []
          for (let j = 0; j < 3; j++)
              vertices.push(new Vertex(vs[f[fm[j]]],f.vertexNormals[j],geom.faceVertexUvs[0][i][j]))
          polys.push(new Polygon(vertices))
      }
  } else if (geom.isBufferGeometry) {
      let vertices, normals, uvs
      let posattr = geom.attributes.position
      let normalattr = geom.attributes.normal
      let uvattr = geom.attributes.uv
      let colorattr = geom.attributes.color
      let index;
      if (geom.index)
          index = geom.index.array;
      else {
          index = new Array((posattr.array.length / posattr.itemSize) | 0);
          for (let i = 0; i < index.length; i++)
              index[i] = i
      }
      let triCount = (index.length / 3) | 0
      polys = new Array(triCount)
      for (let i = 0, pli = 0, l = index.length; i < l; i += 3,
      pli++) {
          let vertices = new Array(3)
          for (let j = 0; j < 3; j++) {
              let vi = index[i + j]
              let vp = vi * 3;
              let vt = vi * 2;
              let x = posattr.array[vp]
              let y = posattr.array[vp + 1]
              let z = posattr.array[vp + 2]
              let nx = normalattr.array[vp]
              let ny = normalattr.array[vp + 1]
              let nz = normalattr.array[vp + 2]
              //let u = uvattr.array[vt]
              //let v = uvattr.array[vt + 1]
              vertices[j] = CSG.Vertex({
                  x,
                  y,
                  z
              },{
                  x: nx,
                  y: ny,
                  z: nz
              });
          }
          polys[pli] = CSG.Polygon(vertices)
      }
  } else
      console.error("Unsupported CSG input type:" + geom.type)
  return CSG.fromPolygons(polys)
}

CSG.toMesh = function (csg, toMatrix, toMaterial) {
  let geom = CSG.toGeometry(csg);
  let inv = new THREE.Matrix4().copy(toMatrix).invert();
  geom.applyMatrix4(inv);
  geom.computeBoundingSphere();
  geom.computeBoundingBox();
  let m = new THREE.Mesh(geom, toMaterial);
  m.matrix.copy(toMatrix);
  m.matrix.decompose(m.position, m.quaternion, m.scale)
  m.rotation.setFromQuaternion(m.quaternion)
  m.updateMatrixWorld();
  m.castShadow = m.receiveShadow = true;
  return m
}

let ttvv0 = new THREE.Vector3()
let tmpm3 = new THREE.Matrix3();
CSG.fromMesh = function(mesh) {
    let csg = CSG.fromGeometry(mesh.geometry)
    tmpm3.getNormalMatrix(mesh.matrix);
    for (let i = 0; i < csg.polygons.length; i++) {
        let p = csg.polygons[i]
        for (let j = 0; j < p.vertices.length; j++) {
            let v = p.vertices[j]
            v.pos.copy(ttvv0.copy(v.pos).applyMatrix4(mesh.matrix));
            v.normal.copy(ttvv0.copy(v.normal).applyMatrix3(tmpm3))
        }
    }
    return csg;
}

// # class Node

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.

CSG.Node = function (polygons) {
  this.plane = null;
  this.front = null;
  this.back = null;
  this.polygons = [];
  if (polygons) this.build(polygons);
};

CSG.Node.prototype = {
  clone: function () {
    var node = new CSG.Node();
    node.plane = this.plane && this.plane.clone();
    node.front = this.front && this.front.clone();
    node.back = this.back && this.back.clone();
    node.polygons = this.polygons.map(function (p) { return p.clone(); });
    return node;
  },

  // Convert solid space to empty space and empty space to solid space.
  invert: function () {
    for (var i = 0; i < this.polygons.length; i++) {
      this.polygons[i].flip();
    }
    this.plane.flip();
    if (this.front) this.front.invert();
    if (this.back) this.back.invert();
    var temp = this.front;
    this.front = this.back;
    this.back = temp;
  },

  // Recursively remove all polygons in `polygons` that are inside this BSP
  // tree.
  clipPolygons: function (polygons) {
    if (!this.plane) return polygons.slice();
    var front = [], back = [];
    for (var i = 0; i < polygons.length; i++) {
      this.plane.splitPolygon(polygons[i], front, back, front, back);
    }
    if (this.front) front = this.front.clipPolygons(front);
    if (this.back) back = this.back.clipPolygons(back);
    else back = [];
    return front.concat(back);
  },

  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `bsp`.
  clipTo: function (bsp) {
    this.polygons = bsp.clipPolygons(this.polygons);
    if (this.front) this.front.clipTo(bsp);
    if (this.back) this.back.clipTo(bsp);
  },

  // Return a list of all polygons in this BSP tree.
  allPolygons: function () {
    var polygons = this.polygons.slice();
    if (this.front) polygons = polygons.concat(this.front.allPolygons());
    if (this.back) polygons = polygons.concat(this.back.allPolygons());
    return polygons;
  },

  // Build a BSP tree out of `polygons`. When called on an existing tree, the
  // new polygons are filtered down to the bottom of the tree and become new
  // nodes there. Each set of polygons is partitioned using the first polygon
  // (no heuristic is used to pick a good split).
  build: function (polygons) {
    if (!polygons.length) return;
    if (!this.plane) this.plane = polygons[0].plane.clone();
    var front = [], back = [];
    for (var i = 0; i < polygons.length; i++) {
      this.plane.splitPolygon(polygons[i], this.polygons, this.polygons, front, back);
    }
    if (front.length) {
      if (!this.front) this.front = new CSG.Node();
      this.front.build(front);
    }
    if (back.length) {
      if (!this.back) this.back = new CSG.Node();
      this.back.build(back);
    }
  }
};
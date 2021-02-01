pico-8 cartridge // http://www.pico-8.com
version 29
__lua__

camera={}
objects={}

function _init()
   camera = perspective2(45, 1, 0.1, 1000)
   add(objects, object(matrix(), cubemesh()))
end

function _update()
end

function _draw()
   cls()
   --drawTriangle(point(2,3), point(8,50), point(30, 6), 12)
   --print(vstr(mapply(rotmatrix(0.25,0,0), vdir(1,0,0))))
   render3d(camera, objects)
end

-------------------------------------------------------------------------------
-- points (vector2d)

function point(x, y)    return {x=x, y=y} end
function padd(p1, p2)   return {x=p1.x+p2.x, y=p1.y+p2.y} end
function psub(p1, p2)   return {x=p1.x-p2.x, y=p1.y-p2.y} end
function pnormsqr(p)    return p.x*p.x + p.y*p.y end
function pdot(p1, p2)   return p1.x * p2.x + p1.y * p2.y end

function sign(p1, p2, p3)
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
end

-- https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
function pointInTriangle(p, t1, t2, t3)
    d1 = sign(p, t1, t2);
    d2 = sign(p, t2, t3);
    d3 = sign(p, t3, t1);

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos);
end

function drawTriangle(p0, p1, p2, col)
   pmin = point(min(min(p0.x, p1.x), p2.x), min(min(p0.y, p1.y), p2.y))
   pmax = point(max(max(p0.x, p1.x), p2.x), max(max(p0.y, p1.y), p2.y))
   for x=pmin.x, pmax.x do
      for y=pmin.y, pmax.y do
         if pointInTriangle(point(x,y), p0, p1, p2) then
            pset(x, y, col)
         end
      end
   end
end

-------------------------------------------------------------------------------
-- vectors

function vpoint(x,y,z) return {x, y, z, 1} end
function vdir(x,y,z) return {x, y, z, 0} end
function vadd(v1, v2) return {v1[1]+v2[1], v1[2]+v2[2], v1[3]+v2[3], v1.w} end
function vsub(v1, v2) return {v1[1]-v2[1], v1[2]-v2[2], v1[3]-v2[3], v1.w} end
function vscale(v1, s) return {v1[1]*s, v1[2]*s, v1[3]*s, v1.w} end
function vdot(v1, v2) return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3] end
function vnormalize(v) return vscale(v, sqrt(vdot(v,v))) end
function vstr(v) return "["..v[1]..", "..v[2]..", "..v[3]..", "..v[4].."]" end

-------------------------------------------------------------------------------
-- matrix

function matrix() return {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}} end
function rotmatrix(p, y, r, m) -- pitch around z, yaw around y, roll around x
   if type(m) == "nil" then
      m = matrix()
   end
   m[1][1] = cos(y)*cos(p)
   m[1][2] = cos(y)*sin(p)*sin(r) - sin(y)*sin(r)
   m[1][3] = cos(y)*sin(p)*cos(r) + sin(y)*sin(r)
   m[2][1] = sin(y)*cos(p)
   m[2][2] = sin(y)*sin(p)*sin(r) + cos(y)*cos(r)
   m[2][3] = sin(y)*sin(p)*cos(r) - cos(y)*sin(r)
   m[3][1] = -sin(p)
   m[3][2] = cos(p)*sin(r)
   m[3][3] = cos(p)*cos(r)
   return m
end
function transmatrix(t, m)
   if type(m) == "nil" then
      m = matrix()
   end
   m[1][4] = t[1]
   m[2][4] = t[2]
   m[3][4] = t[3]
   return m
end
function scalematrix(sx, sy, sz, m)
   if type(m) == "nil" then
      m = matrix()
   end
   m[1][1] *= sx
   m[2][2] *= sy
   m[3][3] *= sz
   m[4][4] = 1
   return m
end
function orth(l, r, t, b, n, f)
   m = matrix()
   m[1][1] = 2 / (r - l)
   m[2][2] = 2 / (t - b)
   m[3][3] = 2 / (f - n)
   m[1][4] = -(r+l) / (r-l)
   m[2][4] = -(t+b) / (t-b)
   m[3][4] = -(f+n) / (f-n)
   m[1][4] = 1
   return m
end
function perspective(l, r, t, b, n, f)
   m = matrix()
   m[1][1] = 2*n / (r - l)
   m[2][2] = 2*n / (t - b)
   m[3][3] = (f+n) / (f-n)
   m[1][3] = -(r+l) / (r-l)
   m[2][3] = -(t+b) / (t-b)
   m[3][4] = -2*f*n / (f-n)
   m[4][3] = 1
   m[4][4] = 0
   return m
end
function perspective2(fovy, aspect, near, far)
   theta = fovy/2
   top = near * sin(theta) / cos(theta) -- tan(x) = sin(x) / cos(x)
   bottom = -top
   right = top * aspect
   left = -right
   return perspective(left, right, bottom, top, near, far);
end
function mmult(m1, m2)
   m = matrix()
   for i=1,4 do
      for j=1,4 do
         m_ij = 0
         for r=1,4 do
            m_ij += m1[i][r] * m2[r][j]
         end
         m[i][j] = m_ij
      end
   end
   return m
end
function mapply(m, v)
   res = vpoint(v[1], v[2], v[3])
   res.w = v.w
   for i=1,4 do
      res[i] = 0
      for j=1,4 do
         res[i] += m[i][j] * v[j]
      end
   end
   return res
end

function meshdata(vertices, textureCoords, normals, triangles)
   return {verts=vertices,
           uvs=textureCoords,
           normals=normals,
           tris=triangles}
end
function object(mat, mesh) return {mat=mat, mesh=mesh} end

function geometryVertexShading(v, uv, n)
   return {v=v, uv=uv, n=n}
end

function geometryProjection(vFrag, objMat, camMat)
   camspaceVert     = mapply(camMat, mapply(objMat, vFrag.v))
   camspaceNormal   = vnormalize(mapply(camMat, mapply(objMat, vFrag.n)))
   return {v=camspaceVert, uv=vFrag.uv, n=camspaceNormal}
end

function geometryClipping(projVertFrag, clippedTriangles)
    -- todo
   add(clippedTriangles, projVertFrag)
end

function geometryScreenMapping(clippedVert, screenMat)
   return mapply(screenMat, clippedVert)
end

function processGeometries(cam, objects)
   screenMat = transmatrix(vpoint(64, 64, 0.5), scalematrix(64, 64, 0.5))
   for obj in all(objects) do
      mesh = obj.mesh
      obj.vFrags = {}
      for t in all(mesh.tris) do
         processedTriangles = {}
         clippedTriangles   = {}
         for tindex=1,3 do
            v = mesh.verts[t[tindex][1]]
            uv = mesh.uvs[t[tindex][2]]
            n = mesh.normals[t[tindex][3]]
            vFrag           = geometryVertexShading(v, uv, n)
            projVertFrag    = geometryProjection(vFrag, obj.mat, cam)
            add(processedTriangles, projVertFrag)
         end
         for processedTriangle in all(processedTriangles) do
            geometryClipping(processedTriangle, clippedTriangles)
         end
         for clippedTriangle in all(clippedTriangles) do
            for tindex=1,3 do
               clippedFrag=clippedTriangles[tindex]
               screenVert = geometryScreenMapping(clippedFrag.v, screenMat)
               clippedFrag.v = screenVert
            end
         end
      end
   end
end
function rasterize()
   -- Rastertrianglesetup
   -- rasterTriangleTraversal
end
function processPixels()
   --pixelShading
   --pixelMerging
end

function render3d(camera, objects)
   processGeometries(camera, objects)
   rasterize()
   processPixels()
end


function cubemesh()
   return meshdata(
      -- vertices
      {vpoint(1.000000, 1.000000, -1.000000),
       vpoint(1.000000, -1.000000, -1.000000),
       vpoint(1.000000, 1.000000, 1.000000),
       vpoint(1.000000, -1.000000, 1.000000),
       vpoint(-1.000000, 1.000000, -1.000000),
       vpoint(-1.000000, -1.000000, -1.000000),
       vpoint(-1.000000, 1.000000, 1.000000),
       vpoint(-1.000000, -1.000000, 1.000000)},
      -- uvs
      {point(0.875000, 0.500000),
       point(0.625000, 0.750000),
       point(0.625000, 0.500000),
       point(0.375000, 1.000000),
       point(0.375000, 0.750000),
       point(0.625000, 0.000000),
       point(0.375000, 0.250000),
       point(0.375000, 0.000000),
       point(0.375000, 0.500000),
       point(0.125000, 0.750000),
       point(0.125000, 0.500000),
       point(0.625000, 0.250000),
       point(0.875000, 0.750000),
       point(0.625000, 1.000000)},
      -- normals
      {vpoint(0.0000, 1.0000, 0.0000),
       vpoint(0.0000, 0.0000, 1.0000),
       vpoint(-1.0000, 0.0000, 0.0000),
       vpoint(0.0000, -1.0000, 0.0000),
       vpoint(1.0000, 0.0000, 0.0000),
       vpoint(0.0000, 0.0000, -1.0000)},
      -- triangles (v,uv,n)
      {{{5,1,1}, {3,2,1}, {1,3,1}},
         {{3,2,2}, {8,4,2}, {4,5,2}},
         {{7,6,3}, {6,7,3}, {8,8,3}},
         {{2,9,4}, {8,10,4}, {6,11,4}},
         {{1,3,5}, {4,5,5}, {2,9,5}},
         {{5,12,6}, {2,9,6}, {6,7,6}},
         {{5,1,1}, {7,13,1}, {3,2,1}},
         {{3,2,2}, {7,14,2}, {8,4,2}},
         {{7,6,3}, {5,12,3}, {6,7,3}},
         {{2,9,4}, {4,5,4}, {8,10,4}},
         {{1,3,5}, {3,2,5}, {4,5,5}},
         {{5,12,6}, {1,3,6}, {2,9,6}}})
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

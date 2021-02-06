pico-8 cartridge // http://www.pico-8.com
version 29
__lua__

-- data
cubeMoveSpeed=0.1

-- runtime globals
camera={}
cube={}
objects={}
cubeRotMat=0
screenMat = 0

-- flow

function _init()
   --camera = perspective2(45, 1, 0.1, 1000)
   --camera = ortho(-5, 5, 5, -5, 0.1, 100)
   camera = perspective(-3, 3, 3, -3, -2, -100)
   cube = object(transmatrix(vpoint(0,0,8)), cubemesh())
   add(objects, cube)

   cubeRotMat= rotmatrix(0.25/30, 0.1/30, 0.05/30)
   screenMat = mmult(transmatrix(vpoint(64, 64, 0.5)), scalematrix(64, 64, 0.5))
end

function _update()
   if (btn(⬆️)) mtranslate(cube.mat, vpoint(0,0,-cubeMoveSpeed))
   if (btn(⬇️)) mtranslate(cube.mat, vpoint(0,0,cubeMoveSpeed))
   for obj in all(objects) do
      local translation = mgettranslation(obj.mat)
      mtranslate(obj.mat, vscale(translation, -1))
      obj.mat = mmult(cubeRotMat, obj.mat)
      mtranslate(obj.mat, translation)
   end
end

function _draw()
   srand(12345)
   cls()
   render3d(camera, objects)
   dflush()
end

-------------------------------------------------------------------------------
-- debug
debugStrs = {}
function dprint(str) add(debugStrs, str) end
function dflush()
   color(7)
   print("mem:"..stat(0).." cpu:"..stat(1).." fps:"..stat(7))
   for s in all(debugStrs) do
      print(s)
   end
   debugStrs = {}
end

-------------------------------------------------------------------------------
-- points (vector2d)

function point2d(x, y)  return {x, y} end
function padd(p1, p2)   return point2d(p1[1]+p2[1], p1[2]+p2[2]) end
function psub(p1, p2)   return point2d(p1[1]-p2[1], p1[2]-p2[2]) end
function pnormsqr(p)    return p[1]*p[1] + p[2]*p[2] end
function pscale(p,s)    return point2d(p[1]*s, p[2]*s) end
function pdot(p1, p2)   return p1[1] * p2[1] + p1[2] * p2[2] end
function pstr(p)        return "<"..p[1]..","..p[2]..">" end

-- https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
function pointInTriangle(p, t1, t2, t3)
   function sign(p1, p2, p3)
      return (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
   end

   local d1 = sign(p, t1, t2);
   local d2 = sign(p, t2, t3);
   local d3 = sign(p, t3, t1);

   local has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
   local has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

   return not (has_neg and has_pos);
end

function drawTriangle(p0, p1, p2, col)
   local pmin = point2d(min(min(p0[1], p1[1]), p2[1]), min(min(p0[2], p1[2]), p2[2]))
   local pmax = point2d(max(max(p0[1], p1[1]), p2[1]), max(max(p0[2], p1[2]), p2[2]))
   for x=pmin[1], pmax[1] do
      for y=pmin[2], pmax[2] do
         if pointInTriangle(point2d(x,y), p0, p1, p2) then
            pset(x, y, col)
         end
      end
   end
end

-------------------------------------------------------------------------------
-- floats

function sign(x)
   if x < 0 then    return -1
   else             return 1
   end
end   
-------------------------------------------------------------------------------
-- vectors

function vpoint(x,y,z) return {x, y, z, 1} end
function vdir(x,y,z) return {x, y, z, 0} end
function vmakedir(v) return {v[1], v[2], v[3], 0} end
function vadd(v1, v2) return {v1[1]+v2[1], v1[2]+v2[2], v1[3]+v2[3], v1[4]+v2[4]} end
function vsub(v1, v2) return {v1[1]-v2[1], v1[2]-v2[2], v1[3]-v2[3], v1.w} end
function vscale(v1, s) return {v1[1]*s, v1[2]*s, v1[3]*s, v1[4]*s} end
function vdot(v1, v2) return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3] end
function vnorm(v) return sqrt(vdot(v,v)) end
function vnormalize(v) local dv = vmakedir(v) return vscale(dv, 1/vnorm(dv)) end
function vdistsqr(v1, v2) local dv = vsub(v2,v1) return vdot(dv,dv) end
function vstr(v) return "["..v[1]..","..v[2]..","..v[3]..","..v[4].."]" end
function vstr3(v) return "["..v[1]..","..v[2]..","..v[3].."]" end
function vcross(v1,v2) return {v1[2]*v2[3] - v1[3]*v2[2],
                               v1[3]*v2[1] - v1[1]*v2[3],
                               v1[1]*v2[2] - v1[2]*v2[1],
                               0}
end

-------------------------------------------------------------------------------
-- matrix

function matrix() return {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}} end
function rotmatrix(x,y,z) -- ammount to rotate around each axis
   function rot_x(xangle)
      local m = matrix()
      m[2][2] = cos(xangle)
      m[2][3] = -sin(xangle)
      m[3][2] = sin(xangle)
      m[3][3] = cos(xangle)
      return m
   end
   function rot_y(yangle)
      local m = matrix()
      m[1][1] = cos(yangle)
      m[1][3] = sin(yangle)
      m[3][1] = -sin(yangle)
      m[3][3] = cos(yangle)
      return m
   end
   function rot_z(zangle)
      local m = matrix()
      m[1][1] = cos(zangle)
      m[1][2] = -sin(zangle)
      m[2][1] = sin(zangle)
      m[2][2] = cos(yangle)
      return m
   end

   return mmult(mmult(rot_x(x), rot_y(y)), rot_z(z))
end
function transmatrix(t)
   local m = matrix()
   m[1][4] = t[1]
   m[2][4] = t[2]
   m[3][4] = t[3]
   return m
end
function scalematrix(sx, sy, sz)
   local m = matrix()
   m[1][1] = sx
   m[2][2] = sy
   m[3][3] = sz
   m[4][4] = 1
   return m
end
function mgettranslation(m)
   return vpoint(m[1][4], m[2][4], m[3][4])
end
function mtranslate(m, t)
   m[1][4] += t[1]
   m[2][4] += t[2]
   m[3][4] += t[3]
end
function ortho(l, r, t, b, n, f)
   local m = matrix()
   m[1][1] = 2 / (r - l)
   m[2][2] = 2 / (t - b)
   m[3][3] = 2 / (f - n)
   m[1][4] = -(r+l) / (r-l)
   m[2][4] = -(t+b) / (t-b)
   m[3][4] = -(f+n) / (f-n)
   m[4][4] = 1
   return m
end
function perspective(l, r, t, b, n, f)
   local m = matrix()
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
   local theta = fovy/2
   local top = near * sin(theta) / cos(theta) -- tan(x) = sin(x) / cos(x)
   local bottom = -top
   local right = top * aspect
   local left = -right
   return perspective(left, right, bottom, top, near, far);
end
function mmult(m1, m2)
   local m = matrix()
   for i=1,4 do
      for j=1,4 do
         local m_ij = 0
         for r=1,4 do
            m_ij += m1[i][r] * m2[r][j]
         end
         m[i][j] = m_ij
      end
   end
   return m
end
function mapply(m, v)
   local res = vpoint(v[1], v[2], v[3])
   res.w = v.w
   for i=1,4 do
      res[i] = 0
      for j=1,4 do
         res[i] += m[i][j] * v[j]
      end
   end
   return res
end
function mprint(m)
   dprint(vstr(m[1]))
   dprint(vstr(m[2]))
   dprint(vstr(m[3]))
   dprint(vstr(m[4]))
end
function meshdata(vertices, textureCoords, normals, triangles, cullBackface)
   return {verts=vertices,
           uvs=textureCoords,
           normals=normals,
           tris=triangles,
           cullBackface=cullBackface}
end
function object(mat, mesh) return {mat=mat, mesh=mesh} end

function fragment(v, uv, n)
   return {vx=v[1], vy=v[2], vz=v[3],
           u=uv[1], v = uv[2],
           nx=n[1], ny=n[2], nz=n[3]}
end
function geometryVertexShading(camMat, objMat, v, uv, n)
   -- projection to camspace
   local camspaceVert   = mapply(camMat, mapply(objMat, v))
   camspaceVert         = vscale(camspaceVert, 1/camspaceVert[4]) // normalize the w component
   local camspaceNormal = vnormalize(mapply(camMat, mapply(objMat, n)))
   return fragment(camspaceVert, uv, camspaceNormal)
end

function geometryClipping(projVertFrag, clippedTriangles)
    -- todo
   add(clippedTriangles, projVertFrag)
end

function geometryScreenMapping(vfrag, screenMat)
   local v = vpoint(vfrag.vx, vfrag.vy, vfrag.vz)
   local screen_v = mapply(screenMat, v)
   --dprint("v: "..vstr(v))
   --dprint("sv: "..vstr(screen_v))
   vfrag.vx = screen_v[1]
   vfrag.vy = screen_v[2]
   vfrag.vz = screen_v[3]
end

function processGeometries(cam, objects)
   -- todo: optimize, only process vertices and normals once
   local vFrags={}
   for obj in all(objects) do
      local mesh = obj.mesh
      for t in all(mesh.tris) do
         local processedTriangleFrags= {}
         local clippedTriangleFrags = {}
         for tindex=1,3 do
            local v     = mesh.verts[t[tindex][1]]
            local uv    = mesh.uvs[t[tindex][2]]
            local n     = mesh.normals[t[tindex][3]]
            local vFrag = geometryVertexShading(cam, obj.mat, v, uv, n)
            --dprint("pv: "..vstr(vFrag.n))
            add(processedTriangleFrags, vFrag)
         end
         for processedTriangle in all(processedTriangleFrags) do
            geometryClipping(processedTriangle, clippedTriangleFrags)
            --dprint(processedTriangle.v[1].." -> "..clippedTriangleFrags[#clippedTriangleFrags].v[1])
         end
         for clippedVertexFrag in all(clippedTriangleFrags) do
            geometryScreenMapping(clippedVertexFrag, screenMat)
            --dprint("frag: "..clippedVertexFrag.vx)
         end
         add(vFrags, clippedTriangleFrags) -- accumulate screen space clipped triangles
      end
      ---------------------------------
      -- local test={}
      -- add(test, vFrags[1])
      -- vFrags = test
      ---------------------------------
   end
   return vFrags
end

function rasterize(vFrags)
   function computeEdgeParams(p1, p2) return {a=-(p2[2]-p1[2]), b=(p2[1]-p1[1]), c=(p2[2]-p1[2])*p1[1] - (p2[1]-p1[1])*p1[2]} end
   function edgeSign(px,py,a,b,c) return a*px + b*py + c end
   function baryInterpVertex(w,u,v,v1,v2,v3) return vadd(vadd(vscale(v1, w), vscale(v2, u)), vscale(v3, v)) end
   function baryInterpPoint(w,u,v,v1,v2,v3) return padd(padd(pscale(v1, w), pscale(v2, u)), pscale(v3, v)) end
   -- rasterTrianglesetup
   -- rasterTriangleTraversal
   local triCount = 0
   local pixelFragCount = 0
   local pFrags = {}
   for i=1,128 do pFrags[i] = {} end

   for vfragTriangle in all(vFrags) do
      triCount += 1
      local tricol = rnd(16)
      local p1 = vpoint(vfragTriangle[1].vx, vfragTriangle[1].vy, vfragTriangle[1].vz)
      local p2 = vpoint(vfragTriangle[2].vx, vfragTriangle[2].vy, vfragTriangle[2].vz)
      local p3 = vpoint(vfragTriangle[3].vx, vfragTriangle[3].vy, vfragTriangle[3].vz)

      -- if obj.mesh.cullBackface then
      --    triNormal = vcross(p1, p2)
      --    if vdot(triNormal, vscale(p1,-1)) < 0 then
      --       break
      --    end
      -- end
      
      local uv1 = point2d(vfragTriangle[1].u, vfragTriangle[1].v)
      local uv2 = point2d(vfragTriangle[2].u, vfragTriangle[2].v)
      local uv3 = point2d(vfragTriangle[3].u, vfragTriangle[3].v)
      local n1 = vpoint(vfragTriangle[1].nx, vfragTriangle[1].ny, vfragTriangle[1].nz)
      local n2 = vpoint(vfragTriangle[2].nx, vfragTriangle[2].ny, vfragTriangle[2].nz)
      local n3 = vpoint(vfragTriangle[3].nx, vfragTriangle[3].ny, vfragTriangle[3].nz)
      local edgeParams = {computeEdgeParams(p1, p2), computeEdgeParams(p2, p3), computeEdgeParams(p3, p1)}
      local pmin_x = flr(min(min(p1[1], p2[1]), p3[1]))
      local pmin_y = flr(min(min(p1[2], p2[2]), p3[2]))
      local pmax_x = flr(max(max(p1[1], p2[1]), p3[1]))
      local pmax_y = flr(max(max(p1[2], p2[2]), p3[2]))
      --dprint("min: "..pmin_x..","..pmin_y.." max: "..pmax_x..","..pmax_y)
      --dprint("params: a:"..edgeParams[1].a.."b: "..edgeParams[1].b.."c: "..edgeParams[1].c)
      for x=pmin_x, pmax_x do
         for y=pmin_y, pmax_y do
            local isInside = true
            local edgeValues = {}
            for e=1,3 do
               -- hack while i figure how to properly do backface culling at the tri level
               edgeValues[e] = edgeSign(x,y, edgeParams[e].a, edgeParams[e].b, edgeParams[e].c)
               --isInside = e==1 or edgeValues[e] == 0 or sign(edgeValues[e]) == sign(edgeValues[e-1])
               isInside = isInside and edgeValues[e] <= 0
               if not isInside then
                  break
               end
            end
            --dprint(edgeValues[1])
            --dprint("("..x..","..y.."): "..edgeValues[1]..","..edgeValues[2]..","..edgeValues[3])
            if isInside then
               local areaSum = (edgeValues[1]  + edgeValues[2] + edgeValues[3])
               -- should use perspective corrected coordinates?
               local barycentric_u = edgeValues[2] / areaSum
               local barycentric_v = edgeValues[3] / areaSum
               local barycentric_w = 1 - barycentric_u - barycentric_v
               -- only interpolate depth instead of whole vertex pos?
               local v = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v, p1, p2, p3)
               local uv = baryInterpPoint(barycentric_w, barycentric_u, barycentric_v, uv1, uv2, uv3)
               local n = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v, n1, n2, n3)
               local pfrag = fragment(v, uv, n)
               pfrag.col = tricol

               -- zbuff check
               if not pFrags[x][y] or pfrag.vz < pFrags[x][y].vz then
                  pFrags[x][y] = pfrag
                  pixelFragCount +=1
               end
            end
         end
      end
   end
   dprint("triangle count: "..triCount)
   dprint("pfrag count: "..pixelFragCount)
   return pFrags
end
function processPixels(pFrags)
   --pixelShading
   --pixelMerging ?
   --pixelRender
   for y=1,128 do
      for x=1,128,8 do
         c1,c2,c3,c4,c5,c6,c7,c8 = 0,0,0,0,0,0,0,0
         if (pFrags[x][y])      c1 = 0xF & flr(pFrags[x][y].col)
         if (pFrags[x+1][y])    c2 = 0xF & flr(pFrags[x+1][y].col)
         if (pFrags[x+2][y])    c3 = 0xF & flr(pFrags[x+2][y].col)
         if (pFrags[x+3][y])    c4 = 0xF & flr(pFrags[x+3][y].col)
         if (pFrags[x+4][y])    c5 = 0xF & flr(pFrags[x+4][y].col)
         if (pFrags[x+5][y])    c6 = 0xF & flr(pFrags[x+5][y].col)
         if (pFrags[x+6][y])    c7 = 0xF & flr(pFrags[x+6][y].col)
         if (pFrags[x+7][y])    c8 = 0xF & flr(pFrags[x+7][y].col)
         -- pset(x-1,y,c1)
         -- pset(x,y,c2)
         -- pset(x+1,y,c3)
         -- pset(x+2,y,c4)
         -- pset(x+3,y,c5)
         -- pset(x+4,y,c6)
         -- pset(x+5,y,c7)
         -- pset(x+6,y,c8)

         local addr = 0x6000+64*(y-1)+x\2
         poke2(addr, (c4 << 12) + (c3 << 8) + (c2 << 4) + c1)
         poke2(addr+2, (c8 << 12) + (c7 << 8) + (c6 << 4) + c5)
         -- local data = (c8<<28)+(c7<<24)+(c6<<20)+(c5<<16)+(c4<<12)+(c3<<8)+(c2<<4)+c1
         -- poke4(addr, data)
      end
   end
end

function render3d(camera, objects)
   local vFrags = processGeometries(camera, objects)
   local pFrags = rasterize(vFrags)
   processPixels(pFrags)
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
      {point2d(0.875000, 0.500000),
       point2d(0.625000, 0.750000),
       point2d(0.625000, 0.500000),
       point2d(0.375000, 1.000000),
       point2d(0.375000, 0.750000),
       point2d(0.625000, 0.000000),
       point2d(0.375000, 0.250000),
       point2d(0.375000, 0.000000),
       point2d(0.375000, 0.500000),
       point2d(0.125000, 0.750000),
       point2d(0.125000, 0.500000),
       point2d(0.625000, 0.250000),
       point2d(0.875000, 0.750000),
       point2d(0.625000, 1.000000)},
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
         {{5,12,6}, {1,3,6}, {2,9,6}}},
   -- cullBackface
   true)
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

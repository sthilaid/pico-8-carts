pico-8 cartridge // http://www.pico-8.com
version 29
__lua__

-------------------------------------------------------------------------------
-- globals

-- data
cubeMoveSpeed=0.1

_pal={8,136,130,133,0,129,1,140,12,11,139,3,131,6,6,5}
red_shades={1,2,3,4,5}
blue_shades={9,8,7,6,5}
green_shades={10,11,12,13,5}
-- color indices for light shading
color_red=0
color_blue=1
color_green=2
bg_color=15

function getShades(c)
   if (c == color_red) return red_shades
   if (c == color_blue) return blue_shades
   if (c == color_green) return green_shades
   return green_shades
end

-- runtime globals
camera=false
cube=false
room=false
l1=false
objects={}
lights={}
cubeRotMat=0
screenMat = 0
campitch=0
camyaw=0
camlen=12

-- runtime debug
triCount = 0
pixelFragCount = 0
updateDT=0
renderStartDT=0
renderEndDT=0
geometryDT=0
rasterDT=0
pixelDT=0

-------------------------------------------------------------------------------
-- flow

function makecam(pitch,yaw, len)
   return mmult(mmult(perspective(-3, 3, 3, -3, -2, -100), transmatrix(vpoint(0,0,len))), rotmatrix(pitch,yaw,0))
end
function _init()
   --camera = perspective2(45, 1, 0.1, 1000)
   --camera = ortho(-5, 5, 5, -5, 0.1, 100)
   camera = makecam(campitch,camyaw,camlen)
   cube = object(mmult(transmatrix(vpoint(0,0,0)), rotmatrix(0.1,0.1,0)), cubemesh())
   room = object(mmult(transmatrix(vpoint(0,0,0)), rotmatrix(0,0,0)), roommesh())
   add(objects, cube)
   add(objects, room)

   l1 = dirlight(vnormalize(vdir(0,-1,-1)), 1)
   add(lights, l1)

   cubeRotMat= rotmatrix(0.25/30, 0.1/30, 0.15/30)
   screenMat = mmult(transmatrix(vpoint(64, 64, 0.5)), scalematrix(64, 64, 0.5))
end

function _update()
   if (btn(⬆️)) camlen -= 0.1
   if (btn(⬇️)) camlen += 0.1
   if (btn(➡️)) camyaw += 0.01
   if (btn(⬅️)) camyaw -= 0.01
   if (btnp(❎)) room.isVisible = not room.isVisible
   camera = makecam(campitch,camyaw,camlen)

   -- local translation = mgettranslation(cube.mat)
   -- mtranslate(cube.mat, vscale(translation, -1))
   -- cube.mat = mmult(cubeRotMat, cube.mat)
   -- mtranslate(cube.mat, translation)

   l1.dir = mapply(cubeRotMat, l1.dir)

   updateDT = stat(1)
end

function _draw()
   srand(12345)
   for i,c in pairs(_pal) do
      pal(i-1,c,1)
   end
   cls(bg_color)
   render3d(camera, objects, lights)
   dflush()
end

-------------------------------------------------------------------------------
-- debug

debugStrs = {}
function dprint(str) add(debugStrs, str) end
function dflush()
   color(14)
   print("mem:"..(stat(0)/2048).."%")
   print("tri: "..triCount.." pfrags: "..pixelFragCount)
   print("g:"..geometryDT.." r:"..rasterDT.." p:"..pixelDT)
   for s in all(debugStrs) do
      print(s)
   end
   debugStrs = {}
   geometryDT, rasterDT, pixelDT, triCount, pixelFragCount = 0,0,0,0,0
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

-- -- https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
-- function pointInTriangle(p, t1, t2, t3)
--    function sign(p1, p2, p3)
--       return (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
--    end

--    local d1 = sign(p, t1, t2);
--    local d2 = sign(p, t2, t3);
--    local d3 = sign(p, t3, t1);

--    local has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
--    local has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

--    return not (has_neg and has_pos);
-- end

-- function drawTriangle(p0, p1, p2, col)
--    local pmin = point2d(min(min(p0[1], p1[1]), p2[1]), min(min(p0[2], p1[2]), p2[2]))
--    local pmax = point2d(max(max(p0[1], p1[1]), p2[1]), max(max(p0[2], p1[2]), p2[2]))
--    for x=pmin[1], pmax[1] do
--       for y=pmin[2], pmax[2] do
--          if pointInTriangle(point2d(x,y), p0, p1, p2) then
--             pset(x, y, col)
--          end
--       end
--    end
-- end

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

-------------------------------------------------------------------------------
-- mesh and objects

function meshdata(vertices, textureCoords, normals, triangles, cullBackface)
   return {verts=vertices,
           uvs=textureCoords,
           normals=normals,
           tris=triangles,
           cullBackface=cullBackface}
end
function object(mat, mesh) return {mat=mat, mesh=mesh, isVisible=true} end

-- lights: 0-directional, 1-point
function dirlight(dir, r0) return {type=0, pos=vpoint(0,0,0), dir=dir, r0=r0, campos=0, camdir=0} end
function pointlight(pos, r0) return {type=1, pos=pos, dir=vdir(0,0,0), r0=r0, campos=0, camdir=0} end

function fragment(v, uv, n, viewdir, color)
   return {vx=v[1], vy=v[2], vz=v[3],
           u=uv[1], v = uv[2],
           nx=n[1], ny=n[2], nz=n[3],
           viewdir=viewdir,
           color=color,
           l1pos=false,
           l1dir=false,
           l1type=false,
           l1r0=false,
   }
end

-------------------------------------------------------------------------------
-- rendering pipeline

function geometryVertexShading(projectionMatrix, mesh, v_index, uv_index, n_index,
                               processedVerticesCache, processedNormalsCache, lights, tcolor)
   local v       = mesh.verts[v_index]

   -- lighting
   local l1dir  = false
   local l1type = false
   local l1r0   = false

   local l1 = lights[1]
   if l1 then
      if (l1.type == 0) l1dir = l1.dir
      if (l1.type == 1) l1dir = vnormalize(vsub(l1.pos, v))
      l1type= l1.type
      l1r0  = l1.r0
   end
   
   -- projection to camspace
   local camspaceVert = processedVerticesCache[v_index]
   if not camspaceVert then
      camspaceVert  = mapply(projectionMatrix, v)
      camspaceVert  = vscale(camspaceVert, 1/camspaceVert[4]) // normalize the w component
   end
   local camspaceNormal  = processedNormalsCache[n_index]
   if not camspaceNormal then
      local n       = mesh.normals[n_index]
      camspaceNormal= vnormalize(mapply(projectionMatrix, n))
   end
   local uv = mesh.uvs[uv_index]

   local viewdir = vscale(vnormalize(vdir(camspaceVert[1],camspaceVert[2],camspaceVert[3])), -1)
   -- fragment creation
   local vfrag = fragment(camspaceVert, uv, camspaceNormal, viewdir, tcolor)
   vfrag.l1dir  = l1dir
   vfrag.l1type = l1type
   vfrag.l1r0   = l1r0
   return vfrag
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

function processGeometries(cam, objects, lights)
   for l in all(lights) do
      l.campos = mapply(cam, l.pos)
      l.camdir = mapply(cam, l.dir)
   end
   for obj in all(objects) do
      if (not obj.isVisible) goto continue
      local mesh = obj.mesh
      local processedVerticesCache = {}
      local processedNormalsCache = {}
      local projectionMatrix = mmult(cam, obj.mat)
      for t in all(mesh.tris) do
         local startDT = stat(1)
         local processedTriangleFrags= {}
         local clippedTriangleFrags = {}
         local triangleColor = t[4]
         for tindex=1,3 do
            local v_index   = t[tindex][1]
            local uv_index  = t[tindex][2]
            local n_index   = t[tindex][3]
            local vFrag = geometryVertexShading(projectionMatrix, mesh, v_index, uv_index, n_index,
                                                processedVerticesCache, processedNormalsCache, lights, triangleColor)
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
         triCount += 1
         geometryDT += stat(1) - startDT
         rasterizeTriangle(clippedTriangleFrags)
         --if (triCount == 8) return
      end
      ::continue::
   end
end

function computeEdgeParams(p1, p2) return {a=-(p2[2]-p1[2]), b=(p2[1]-p1[1]), c=(p2[2]-p1[2])*p1[1] - (p2[1]-p1[1])*p1[2]} end
function edgeSign(px,py,a,b,c) return a*px + b*py + c end
function baryInterpVertex(w,u,v,v1,v2,v3) return vadd(vadd(vscale(v1, w), vscale(v2, u)), vscale(v3, v)) end
function baryInterpPoint(w,u,v,v1,v2,v3) return padd(padd(pscale(v1, w), pscale(v2, u)), pscale(v3, v)) end
function baryInterpValue(w,u,v,v1,v2,v3) return v1*w + v2*u + v3*v end
function setPixelFragment(tbl, x, y, frag) tbl[(128*y)+x] = frag end
function getPixelFragment(tbl, x, y) return tbl[(128*y)+x] end

function rasterizeTriangle(clippedTriangleFrags)
   local tricol = clippedTriangleFrags[1].color --rnd(16)
   local p1 = vpoint(clippedTriangleFrags[1].vx, clippedTriangleFrags[1].vy, clippedTriangleFrags[1].vz)
   local p2 = vpoint(clippedTriangleFrags[2].vx, clippedTriangleFrags[2].vy, clippedTriangleFrags[2].vz)
   local p3 = vpoint(clippedTriangleFrags[3].vx, clippedTriangleFrags[3].vy, clippedTriangleFrags[3].vz)

   -- if obj.mesh.cullBackface then
   --    triNormal = vcross(p1, p2)
   --    if vdot(triNormal, vscale(p1,-1)) < 0 then
   --       break
   --    end
   -- end
   
   local uv1 = point2d(clippedTriangleFrags[1].u, clippedTriangleFrags[1].v)
   local uv2 = point2d(clippedTriangleFrags[2].u, clippedTriangleFrags[2].v)
   local uv3 = point2d(clippedTriangleFrags[3].u, clippedTriangleFrags[3].v)
   local n1 = vpoint(clippedTriangleFrags[1].nx, clippedTriangleFrags[1].ny, clippedTriangleFrags[1].nz)
   local n2 = vpoint(clippedTriangleFrags[2].nx, clippedTriangleFrags[2].ny, clippedTriangleFrags[2].nz)
   local n3 = vpoint(clippedTriangleFrags[3].nx, clippedTriangleFrags[3].ny, clippedTriangleFrags[3].nz)
   local edgeParams = {computeEdgeParams(p1, p2), computeEdgeParams(p2, p3), computeEdgeParams(p3, p1)}
   local pmin_x = flr(min(min(p1[1], p2[1]), p3[1]))
   local pmin_y = flr(min(min(p1[2], p2[2]), p3[2]))
   local pmax_x = flr(max(max(p1[1], p2[1]), p3[1]))
   local pmax_y = flr(max(max(p1[2], p2[2]), p3[2]))
   --dprint("min: "..pmin_x..","..pmin_y.." max: "..pmax_x..","..pmax_y)
   --dprint("params: a:"..edgeParams[1].a.."b: "..edgeParams[1].b.."c: "..edgeParams[1].c)
   for x=pmin_x, pmax_x do
      for y=pmin_y, pmax_y do
         local pixelStartDT = stat(1)
         local isInside = true
         local edgeValues = {}
         for e=1,3 do
            edgeValues[e] = edgeSign(x,y, edgeParams[e].a, edgeParams[e].b, edgeParams[e].c)
            --isInside = e==1 or edgeValues[e] == 0 or sign(edgeValues[e]) == sign(edgeValues[e-1])
            -- hack while i figure how to properly do backface culling at the tri level
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
            local depth  = baryInterpValue(barycentric_w, barycentric_u, barycentric_v, p1[3], p2[3], p3[3])
            local uv     = baryInterpPoint(barycentric_w, barycentric_u, barycentric_v, uv1, uv2, uv3)
            local n      = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v, n1, n2, n3)
            local l1dir  = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v,
                                            clippedTriangleFrags[1].l1dir,
                                            clippedTriangleFrags[2].l1dir,
                                            clippedTriangleFrags[3].l1dir)
            local l1type = clippedTriangleFrags[1].l1type
            local l1r0   = clippedTriangleFrags[1].l1r0
            local viewdir= baryInterpVertex(barycentric_w, barycentric_u, barycentric_v,
                                            clippedTriangleFrags[1].viewdir,
                                            clippedTriangleFrags[2].viewdir,
                                            clippedTriangleFrags[3].viewdir)
            local pfrag  = fragment(vpoint(x,y,depth), uv, n, viewdir, tricol)
            pfrag.l1dir  = l1dir
            pfrag.l1type = l1type
            pfrag.l1r0   = l1r0
            rasterDT += stat(1) - pixelStartDT
            processPixelFragment(pfrag)
         end
      end
   end
end

function processPixelFragment(pfrag)
   if (pget(pfrag.vx, pfrag.vy) != bg_color) return // kind of depth test
   
   local startDT = stat(1)
   pixelFragCount += 1

   -- lighting
   local lightRatio = max(0,vdot(pfrag.l1dir, vdir(pfrag.nx,pfrag.ny, pfrag.nz)))
   local color_shades = getShades(pfrag.color)
   local shadedColorIndex = min(flr(lightRatio * (#color_shades-1))+1, #color_shades)
   local shadedColor = color_shades[shadedColorIndex]
   
   -- pixel render
   pset(pfrag.vx, pfrag.vy, shadedColor)
   pixelDT += stat(1) - startDT
end

function render3d(camera, objects, lights)
   renderStartDT = stat(1)
   processGeometries(camera, objects, lights)
   renderEndDT = stat(1)
end

-------------------------------------------------------------------------------
-- cube mesh data (from blender .obj export)

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
      {vdir(0.0000, 1.0000, 0.0000),
       vdir(0.0000, 0.0000, 1.0000),
       vdir(-1.0000, 0.0000, 0.0000),
       vdir(0.0000, -1.0000, 0.0000),
       vdir(1.0000, 0.0000, 0.0000),
       vdir(0.0000, 0.0000, -1.0000)},
      -- triangles (v,uv,n) - col
      {{{5,1,1}, {3,2,1}, {1,3,1},      color_red},
         {{3,2,2}, {8,4,2}, {4,5,2},    color_red},
         {{7,6,3}, {6,7,3}, {8,8,3},    color_red},
         {{2,9,4}, {8,10,4}, {6,11,4},  color_red},
         {{1,3,5}, {4,5,5}, {2,9,5},    color_red},
         {{5,12,6}, {2,9,6}, {6,7,6},   color_red},
         {{5,1,1}, {7,13,1}, {3,2,1},   color_red},
         {{3,2,2}, {7,14,2}, {8,4,2},   color_red},
         {{7,6,3}, {5,12,3}, {6,7,3},   color_red},
         {{2,9,4}, {4,5,4}, {8,10,4},   color_red},
         {{1,3,5}, {3,2,5}, {4,5,5},    color_red},
         {{5,12,6}, {1,3,6}, {2,9,6},   color_red}},
   -- cullBackface
   true)
end
function roommesh()
   return meshdata(
        {vpoint(5.000000, 5.000000, -5.000000),
         vpoint(5.000000, -5.000000, -5.000000),
         vpoint(5.000000, 5.000000, 5.000000),
         vpoint(5.000000, -5.000000, 5.000000),
         vpoint(-5.000000, 5.000000, -5.000000),
         vpoint(-5.000000, -5.000000, -5.000000),
         vpoint(-5.000000, 5.000000, 5.000000),
         vpoint(-5.000000, -5.000000, 5.000000)},
        {point2d(0.375000, 0.750000),
         point2d(0.375000, 1.000000),
         point2d(0.625000, 1.000000),
         point2d(0.625000, 0.750000),
         point2d(0.375000, 0.000000),
         point2d(0.375000, 0.250000),
         point2d(0.625000, 0.250000),
         point2d(0.625000, 0.000000),
         point2d(0.125000, 0.500000),
         point2d(0.125000, 0.750000),
         point2d(0.375000, 0.500000),
         point2d(0.625000, 0.500000)},
        {vdir(0.0000, 0.0000, -1.0000),
         vdir(1.0000, 0.0000, 0.0000),
         vdir(0.0000, 1.0000, 0.0000),
         vdir(-1.0000, 0.0000, 0.0000)},
        {{{8,1,1}, {3,2,1}, {4,3,1},    color_blue},
         {{6,4,2}, {7,5,2}, {8,6,2},    color_blue},
         {{8,7,3}, {2,8,3}, {6,9,3},    color_blue},
         {{4,3,4}, {1,10,4}, {2,8,4},   color_blue},
         {{8,1,1}, {7,11,1}, {3,2,1},   color_blue},
         {{6,4,2}, {5,12,2}, {7,5,2},   color_blue},
         {{8,7,3}, {4,3,3}, {2,8,3},    color_blue},
         {{4,3,4}, {3,2,4}, {1,10,4},   color_blue}}
        )
end


__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

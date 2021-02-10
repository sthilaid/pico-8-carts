pico-8 cartridge // http://www.pico-8.com
version 29
__lua__

-------------------------------------------------------------------------------
-- globals

-- data
cubeMoveSpeed=0.1

_pal={0,128,133,5,134,6,135,15,7}
--_pal={8,136,132,128,0,129,1,140,12,11,139,3,131,5,134,6}
-- red_shades={4,3,2,1,0}
-- blue_shades={4,6,6,7,8}
-- green_shades={4,12,11,10,9}
-- white_shades={4,13,14,15}
white_shades={0,1,2,3,4,5,6,7,8}
-- color indices for light shading
color_red=0
color_blue=1
color_green=2
color_white=3
bg_color=11

function getShades(c)
   if (c == color_red) return red_shades
   if (c == color_blue) return blue_shades
   if (c == color_green) return green_shades
   return white_shades
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
camlen=5
shouldRotateLight=true

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
   rabbit = object(matrix(), rabbitMesh())
   -- add(objects, cube)
   -- add(objects, room)
   add(objects, rabbit)

   l1 = dirlight(vnormalize(vdir(0,0,-1)), 0.0)
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
   if (btnp(🅾️)) shouldRotateLight = not shouldRotateLight
   camera = makecam(campitch,camyaw,camlen)

   -- local translation = mgettranslation(cube.mat)
   -- mtranslate(cube.mat, vscale(translation, -1))
   -- cube.mat = mmult(cubeRotMat, cube.mat)
   -- mtranslate(cube.mat, translation)

   if (shouldRotateLight) l1.dir = mapply(rotmatrix(0,0.01,0), l1.dir)

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
   color(15)
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
function clamp(min,max,x)
   if x < min then return min
   elseif x > max then return max
   else return x
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
   local v      = mesh.verts[v_index]

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
   local lightRatio = 1.0
   if pfrag.l1type == 0 then
      lightRatio = clamp(0,1,pfrag.l1r0+max(0,vdot(pfrag.l1dir, vdir(pfrag.nx,pfrag.ny, pfrag.nz))))
   end
   local color_shades = getShades(pfrag.color)
   local shadedColorIndex = min(flr(lightRatio * (#color_shades-1))+1, #color_shades)
   local shadedColor = color_shades[shadedColorIndex]
   --print(shadedColor)
   
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

rabbitVertices = {
   vpoint(0.000000, -1.391679, -0.797764),
   vpoint(-0.135392, 0.611448, -0.507234),
   vpoint(0.517557, 0.609942, 0.052065),
   vpoint(0.056643, 0.604025, 0.500860),
   vpoint(-0.473639, 0.610471, 0.213351),
   vpoint(0.398882, -1.391679, -0.690884),
   vpoint(0.690884, -1.391679, -0.398882),
   vpoint(0.797764, -1.391679, 0.000000),
   vpoint(0.690884, -1.391679, 0.398882),
   vpoint(0.398882, -1.391679, 0.690884),
   vpoint(0.000000, -1.391679, 0.797764),
   vpoint(-0.398882, -1.391679, 0.690884),
   vpoint(-0.690883, -1.391679, 0.398882),
   vpoint(-0.797764, -1.391679, 0.000000),
   vpoint(-0.690884, -1.391679, -0.398881),
   vpoint(-0.398882, -1.391679, -0.690883),
   vpoint(0.007565, 1.280805, -0.534633),
   vpoint(-0.001319, 0.824538, -0.592229),
   vpoint(0.310472, 0.972387, -0.526726),
   vpoint(-0.004148, 0.435775, -0.283392),
   vpoint(-0.085475, 1.602911, -0.022652),
   vpoint(0.457809, 1.236295, -0.356192),
   vpoint(0.386913, 0.649667, -0.355688),
   vpoint(0.149670, 1.545669, -0.219466),
   vpoint(0.619546, 1.066289, 0.001204),
   vpoint(0.504853, 0.650576, 0.163492),
   vpoint(0.517314, 1.077296, 0.298671),
   vpoint(0.330113, 1.505936, 0.119562),
   vpoint(0.300922, 0.995839, 0.532239),
   vpoint(0.246869, 0.529812, 0.342540),
   vpoint(0.157971, 0.390043, -0.016812),
   vpoint(0.129921, 1.364891, 0.469963),
   vpoint(-0.000000, 0.963595, 0.612501),
   vpoint(-0.133589, 0.587250, 0.454342),
   vpoint(-0.162013, 1.455557, 0.363314),
   vpoint(-0.306251, 1.004631, 0.530442),
   vpoint(-0.117633, 0.423134, 0.210932),
   vpoint(-0.519327, 1.063339, 0.300116),
   vpoint(-0.402400, 1.433085, 0.146065),
   vpoint(-0.615259, 1.078989, 0.007206),
   vpoint(-0.497713, 0.646425, 0.157127),
   vpoint(-0.326746, 1.373354, -0.349210),
   vpoint(-0.537144, 0.901938, -0.308731),
   vpoint(-0.267728, 0.433814, -0.075064),
   vpoint(-0.299692, 1.063339, -0.519572),
   vpoint(-0.273868, 0.687422, -0.459224),
   vpoint(0.403289, 2.117685, 0.035751),
   vpoint(-0.242199, 2.117975, -0.073327),
   vpoint(0.212305, 2.118290, -0.061642),
   vpoint(-0.357801, 2.117975, 0.073328),
   vpoint(-0.368098, 1.280821, -0.110259),
   vpoint(-0.102704, 1.264386, -0.020220),
   vpoint(-0.095093, 1.641306, -0.080582),
   vpoint(-0.045583, 1.810991, 0.025934),
   vpoint(-0.231820, 1.283460, 0.110408),
   vpoint(-0.504907, 1.641307, 0.080582),
   vpoint(-0.497220, 1.264120, 0.020212),
   vpoint(-0.554417, 1.810991, -0.025934),
   vpoint(0.231902, 1.280821, -0.110259),
   vpoint(0.532326, 1.829810, -0.061515),
   vpoint(0.497296, 1.264386, -0.020220),
   vpoint(0.552316, 1.641104, 0.024805),
   vpoint(0.354862, 1.266797, 0.089824),
   vpoint(0.422333, 1.878187, 0.116543),
   vpoint(0.073201, 1.844584, 0.066144),
   vpoint(0.102780, 1.264120, 0.020212),
   vpoint(0.047684, 1.641104, -0.024805),
   vpoint(0.007158, 0.302965, -0.779088),
   vpoint(0.388830, 0.381708, -0.651643),
   vpoint(0.666557, 0.363844, -0.373362),
   vpoint(0.777647, 0.311073, 0.010863),
   vpoint(0.668031, 0.311074, 0.398231),
   vpoint(0.376527, 0.374131, 0.661485),
   vpoint(-0.010862, 0.311073, 0.777647),
   vpoint(-0.390574, 0.363844, 0.656619),
   vpoint(-0.678894, 0.311074, 0.379417),
   vpoint(-0.777647, 0.311074, -0.010862),
   vpoint(-0.656824, 0.374131, -0.384599),
   vpoint(-0.379417, 0.311074, -0.678893),
   vpoint(0.112758, 1.095153, 0.473684),
   vpoint(0.189389, 1.097969, 0.375961),
   vpoint(0.231078, 1.039332, 0.486629),
   vpoint(0.212252, 1.143772, 0.593166),
   vpoint(0.321897, 1.145405, 0.438386),
   vpoint(0.253915, 1.236605, 0.519741),
   vpoint(0.171172, 1.241927, 0.423361),
   vpoint(-0.212785, 1.132481, 0.594174),
   vpoint(-0.311047, 1.144340, 0.448341),
   vpoint(-0.202055, 1.165019, 0.358402),
   vpoint(-0.218891, 1.251477, 0.486628),
   vpoint(-0.201075, 1.033039, 0.473684),
   vpoint(-0.085836, 1.169609, 0.451970)}

rabbitUVs = {
    point2d(0.529308, 0.976941),
    point2d(0.328208, 0.969205),
    point2d(0.416667, 0.846180),
    point2d(0.750000, 0.661596),
    point2d(0.604218, 0.627158),
    point2d(0.666667, 0.494235),
    point2d(0.750000, 0.416675),
    point2d(0.615460, 0.318621),
    point2d(0.751678, 0.158505),
    point2d(0.578685, 0.075052),
    point2d(0.500000, 0.541986),
    point2d(0.450202, 0.321009),
    point2d(0.351973, 0.237748),
    point2d(0.204358, 0.289258),
    point2d(0.177477, 0.143041),
    point2d(0.488493, 0.095560),
    point2d(0.333333, 0.505851),
    point2d(0.291173, 0.707252),
    point2d(0.250000, 0.489589),
    point2d(0.120417, 0.926460),
    point2d(0.184578, 0.773072),
    point2d(0.454815, 0.784807),
    point2d(0.125000, 0.000000),
    point2d(0.030782, 0.114817),
    point2d(0.053751, 0.752787),
    point2d(0.166667, 0.510411),
    point2d(0.083333, 0.541803),
    point2d(0.050310, 0.318860),
    point2d(1.000000, 0.797371),
    point2d(0.846001, 0.938715),
    point2d(0.872447, 0.708921),
    point2d(1.000000, 0.548867),
    point2d(0.833333, 0.541804),
    point2d(0.916667, 0.458153),
    point2d(0.875000, 0.000000),
    point2d(0.955693, 0.153044),
    point2d(0.833333, 0.338405),
    point2d(0.083333, 0.500000),
    point2d(0.001673, 0.825168),
    point2d(0.000000, 0.500000),
    point2d(0.951713, 0.996829),
    point2d(0.694120, 1.000000),
    point2d(0.805398, 0.937337),
    point2d(0.721750, 0.931242),
    point2d(0.861967, 0.875090),
    point2d(0.848491, 0.997802),
    point2d(0.549573, 0.756544),
    point2d(0.766741, 0.750806),
    point2d(0.221750, 0.931242),
    point2d(0.361967, 0.875090),
    point2d(0.394304, 0.997537),
    point2d(0.049643, 0.755641),
    point2d(0.266802, 0.750807),
    point2d(0.582861, 0.941658),
    point2d(0.321715, 0.940209),
    point2d(0.186407, 0.995952),
    point2d(0.049643, 0.755641),
    point2d(0.218930, 0.875091),
    point2d(0.266802, 0.750807),
    point2d(0.166667, 0.851921),
    point2d(0.250000, 0.826127),
    point2d(0.500000, 0.826127),
    point2d(0.583333, 0.849138),
    point2d(0.666667, 0.826127),
    point2d(0.745536, 0.969766),
    point2d(0.974065, 0.970648),
    point2d(0.916667, 0.855045),
    point2d(0.998327, 0.825168),
    point2d(0.867676, 0.225689),
    point2d(0.636363, 0.078731),
    point2d(0.589503, 0.231850),
    point2d(0.272727, 0.000000),
    point2d(0.227571, 0.116845),
    point2d(0.379115, 0.131146),
    point2d(0.560647, 0.446067),
    point2d(0.304455, 0.348877),
    point2d(0.954545, 0.393651),
    point2d(0.841391, 0.386401),
    point2d(0.173071, 0.168881),
    point2d(0.027841, 0.205683),
    point2d(0.390374, 0.174953),
    point2d(0.477994, 0.040616),
    point2d(0.067873, 0.039899),
    point2d(0.868987, 0.226450),
    point2d(0.636363, 0.078731),
    point2d(0.618942, 0.265918),
    point2d(0.750307, 0.432483),
    point2d(0.590610, 0.355537),
    point2d(0.411566, 0.266699),
    point2d(0.545454, 0.472382),
    point2d(0.211872, 0.238265),
    point2d(0.340187, 0.431766),
    point2d(0.066589, 0.248231),
    point2d(0.181818, 0.393651),
    point2d(0.455921, 0.052826),
    point2d(0.272727, 0.078731),
    point2d(0.089815, 0.052339),
    point2d(0.916667, 0.500000),
    point2d(0.833333, 0.846180),
    point2d(0.833333, 0.500000),
    point2d(0.750000, 0.826127),
    point2d(0.750000, 0.500000),
    point2d(0.666667, 0.500000),
    point2d(0.583333, 0.500000),
    point2d(0.500000, 0.500000),
    point2d(0.416667, 0.500000),
    point2d(0.333333, 0.826127),
    point2d(0.333333, 0.500000),
    point2d(0.250000, 0.500000),
    point2d(0.166667, 0.500000),
    point2d(0.310000, 0.026077),
    point2d(0.206885, 0.449368),
    point2d(0.080855, 0.115455),
    point2d(0.083333, 0.826127),
    point2d(0.750000, 0.010000),
    point2d(0.542154, 0.370000),
    point2d(0.957846, 0.370000),
    point2d(0.416667, 0.549337),
    point2d(0.052500, 0.976740),
    point2d(0.655012, 0.858998),
    point2d(0.000000, 0.548867),
    point2d(1.000000, 0.298899),
    point2d(0.084337, 0.999216),
    point2d(0.549884, 0.750166),
    point2d(0.999961, 0.755631),
    point2d(0.718930, 0.875091),
    point2d(0.766741, 0.750806),
    point2d(0.999961, 0.755631),
    point2d(0.463715, 0.282519),
    point2d(1.000000, 0.500000),
    point2d(0.990000, 0.250000),
    point2d(0.957846, 0.130000),
    point2d(0.870000, 0.042154),
    point2d(0.630000, 0.042154),
    point2d(0.542154, 0.130000),
    point2d(0.510000, 0.250000),
    point2d(0.630000, 0.457846),
    point2d(0.750000, 0.490000),
    point2d(0.870000, 0.457846)}

rabbitNormals = {
    vpoint(-0.2645, 0.8208, 0.5063),
    vpoint(0.3750, 0.3463, -0.8599),
    vpoint(0.1924, -0.6113, -0.7677),
    vpoint(0.1488, 0.1210, -0.9814),
    vpoint(0.3887, -0.8392, -0.3803),
    vpoint(0.8856, -0.1074, -0.4519),
    vpoint(0.5229, -0.8301, 0.1935),
    vpoint(0.0398, -0.8273, 0.5604),
    vpoint(0.1560, -0.9330, 0.3243),
    vpoint(0.2225, 0.2616, 0.9392),
    vpoint(0.0889, 0.9245, 0.3706),
    vpoint(-0.1047, -0.9943, 0.0178),
    vpoint(-0.3349, 0.9005, 0.2774),
    vpoint(-0.7416, -0.2517, 0.6219),
    vpoint(-0.5168, -0.8217, 0.2405),
    vpoint(-0.5378, 0.8233, -0.1814),
    vpoint(-0.8065, 0.5596, -0.1907),
    vpoint(-0.7362, 0.2755, -0.6182),
    vpoint(-0.0856, -0.9894, -0.1177),
    vpoint(-0.7051, -0.5969, -0.3828),
    vpoint(-0.3608, -0.1719, -0.9166),
    vpoint(-0.2588, 0.0117, -0.9659),
    vpoint(-0.3445, 0.4289, -0.8351),
    vpoint(0.4067, 0.4508, -0.7946),
    vpoint(0.6958, 0.2251, -0.6820),
    vpoint(0.6952, -0.1323, 0.7065),
    vpoint(-0.6958, 0.2251, 0.6820),
    vpoint(-0.6960, -0.1323, -0.7057),
    vpoint(0.8551, 0.2255, 0.4669),
    vpoint(-0.1647, 0.3035, 0.9385),
    vpoint(-0.3443, 0.0464, -0.9377),
    vpoint(-0.6867, -0.1843, -0.7032),
    vpoint(-0.6645, 0.7419, -0.0897),
    vpoint(-0.1282, 0.6966, 0.7060),
    vpoint(0.5306, 0.7419, 0.4099),
    vpoint(0.1004, 0.6835, -0.7230),
    vpoint(0.7492, -0.4020, 0.5265),
    vpoint(-0.3755, -0.8700, -0.3195),
    vpoint(0.5212, 0.7507, -0.4060),
    vpoint(-0.5552, 0.6548, 0.5128),
    vpoint(-0.7656, 0.4592, 0.4506),
    vpoint(0.8121, 0.0913, 0.5763),
    vpoint(0.3019, 0.3324, -0.8935),
    vpoint(-0.7842, 0.1023, -0.6120),
    vpoint(0.4982, -0.6744, -0.5450),
    vpoint(-0.3923, -0.6906, 0.6077),
    vpoint(0.7021, -0.5142, 0.4926),
    vpoint(0.5221, 0.5864, 0.6192),
    vpoint(0.3616, 0.7943, -0.4881),
    vpoint(-0.5500, 0.6572, -0.5154),
    vpoint(-0.7284, 0.4384, 0.5265),
    vpoint(0.5355, -0.5533, -0.6380),
    vpoint(-0.4435, -0.5915, -0.6734),
    vpoint(-0.6919, -0.5884, 0.4184),
    vpoint(0.7070, 0.0201, -0.7070),
    vpoint(0.9658, 0.0131, -0.2588),
    vpoint(0.9621, 0.0130, 0.2723),
    vpoint(0.6726, 0.0208, 0.7397),
    vpoint(0.2588, 0.0131, 0.9658),
    vpoint(-0.3011, 0.0200, 0.9534),
    vpoint(-0.7070, 0.0131, 0.7070),
    vpoint(-0.9694, 0.0130, 0.2453),
    vpoint(-0.9523, 0.0208, -0.3044),
    vpoint(0.0091, 0.9999, 0.0056),
    vpoint(-0.7070, 0.0131, -0.7070),
    vpoint(0.3124, 0.0228, -0.9497),
    vpoint(0.0000, -1.0000, -0.0000),
    vpoint(0.2961, 0.8914, -0.3433),
    vpoint(0.6368, 0.5171, 0.5719),
    vpoint(-0.4453, 0.8709, -0.2078),
    vpoint(0.3880, 0.8258, 0.4093),
    vpoint(0.3368, 0.6405, -0.6901),
    vpoint(0.3559, -0.3704, -0.8580),
    vpoint(0.6291, 0.7365, -0.2485),
    vpoint(0.8243, -0.1001, -0.5573),
    vpoint(0.6411, -0.7537, -0.1443),
    vpoint(0.8183, 0.5659, -0.1011),
    vpoint(0.9190, -0.3348, -0.2082),
    vpoint(0.2299, 0.9732, -0.0083),
    vpoint(0.9361, -0.1308, 0.3265),
    vpoint(0.8361, 0.4778, 0.2697),
    vpoint(0.6393, -0.3524, 0.6835),
    vpoint(0.2761, -0.3897, 0.8786),
    vpoint(0.7527, -0.2187, 0.6210),
    vpoint(0.1211, 0.8950, 0.4294),
    vpoint(-0.2037, 0.3969, 0.8950),
    vpoint(-0.2835, -0.2843, 0.9159),
    vpoint(-0.6184, 0.4608, 0.6366),
    vpoint(-0.6209, -0.3817, 0.6847),
    vpoint(-0.1839, 0.8629, -0.4707),
    vpoint(-0.4969, -0.7343, 0.4625),
    vpoint(-0.9417, -0.1518, 0.3003),
    vpoint(-0.8644, 0.4002, 0.3045),
    vpoint(-0.8483, 0.3440, -0.4025),
    vpoint(-0.7808, -0.5733, -0.2484),
    vpoint(-0.9552, -0.2859, -0.0760),
    vpoint(-0.3915, -0.7708, -0.5026),
    vpoint(-0.5885, -0.1674, -0.7910),
    vpoint(-0.1378, 0.1267, -0.9823),
    vpoint(0.0867, -0.9954, -0.0416),
    vpoint(-0.0693, -0.6202, -0.7814),
    vpoint(0.6424, 0.5752, 0.5064),
    vpoint(-0.6424, 0.5752, -0.5064),
    vpoint(0.0359, -0.9951, -0.0918),
    vpoint(0.0628, 0.0346, -0.9974),
    vpoint(0.3086, -0.1565, -0.9382),
    vpoint(0.9430, -0.0711, -0.3251),
    vpoint(0.2231, 0.0769, 0.9718),
    vpoint(-0.0634, 0.0347, 0.9974),
    vpoint(-0.3076, -0.1565, 0.9386),
    vpoint(-0.9431, -0.0712, 0.3249),
    vpoint(-0.2226, 0.0764, -0.9719),
    vpoint(-0.3802, 0.5428, 0.7489),
    vpoint(0.5994, -0.1807, 0.7798),
    vpoint(0.0537, 0.0591, -0.9968),
    vpoint(0.3151, -0.0885, -0.9449),
    vpoint(0.9130, -0.0857, -0.3988),
    vpoint(0.8470, 0.2927, 0.4437),
    vpoint(-0.2642, -0.0895, 0.9603),
    vpoint(-0.8988, -0.0799, 0.4309),
    vpoint(-0.8862, 0.2771, -0.3713),
    vpoint(-0.1892, 0.5986, -0.7784),
    vpoint(-0.4477, 0.6892, -0.5697),
    vpoint(-0.7568, 0.6250, 0.1915),
    vpoint(-0.5792, 0.6614, 0.4765),
    vpoint(0.1115, 0.6688, 0.7350),
    vpoint(0.7514, 0.6246, 0.2126),
    vpoint(0.7401, 0.6611, -0.1232),
    vpoint(-0.2507, 0.0097, -0.9680),
    vpoint(0.7083, 0.0196, -0.7056),
    vpoint(0.9611, 0.0173, -0.2755),
    vpoint(0.9659, 0.0098, 0.2588),
    vpoint(0.7071, 0.0098, 0.7071),
    vpoint(0.2842, 0.0196, 0.9586),
    vpoint(-0.2588, 0.0098, 0.9659),
    vpoint(-0.6946, 0.0173, 0.7192),
    vpoint(-0.9659, 0.0098, 0.2588),
    vpoint(-0.9659, 0.0098, -0.2588),
    vpoint(-0.0048, 1.0000, 0.0083),
    vpoint(-0.7254, 0.0196, -0.6880),
    vpoint(0.2588, 0.0096, -0.9659),
    vpoint(0.4034, 0.8462, -0.3482),
    vpoint(0.6008, 0.3956, 0.6946),
    vpoint(0.1921, -0.4374, 0.8785),
    vpoint(-0.1909, 0.3829, 0.9039),
    vpoint(-0.6007, 0.4396, 0.6678),
    vpoint(-0.1190, 0.7859, -0.6068),
    vpoint(-0.0830, -0.9945, 0.0632),
    vpoint(-0.0481, -0.9980, -0.0404),
    vpoint(0.4603, -0.0894, 0.8832),
    vpoint(-0.1402, -0.0278, 0.9897)}

rabbitTriangles = {
    {{4,1,1}, {5,2,1}, {75,3,1}, color_white},
    {{17,4,2}, {22,5,2}, {19,6,2}, color_white},
    {{18,7,3}, {23,8,3}, {20,9,3}, color_white},
    {{17,4,4}, {19,6,4}, {18,7,4}, color_white},
    {{20,9,5}, {23,8,5}, {31,10,5}, color_white},
    {{22,5,6}, {25,11,6}, {23,8,6}, color_white},
    {{26,12,7}, {30,13,7}, {31,10,7}, color_white},
    {{30,13,8}, {34,14,8}, {37,15,8}, color_white},
    {{30,13,9}, {37,15,9}, {31,16,9}, color_white},
    {{29,17,10}, {32,18,10}, {33,19,10}, color_white},
    {{21,20,11}, {35,21,11}, {28,22,11}, color_white},
    {{31,23,12}, {37,15,12}, {44,24,12}, color_white},
    {{35,21,13}, {21,20,13}, {39,25,13}, color_white},
    {{36,26,14}, {38,27,14}, {41,28,14}, color_white},
    {{41,28,15}, {44,24,15}, {37,15,15}, color_white},
    {{39,29,16}, {21,30,16}, {42,31,16}, color_white},
    {{39,29,17}, {42,31,17}, {40,32,17}, color_white},
    {{42,31,18}, {45,33,18}, {43,34,18}, color_white},
    {{31,35,19}, {44,36,19}, {20,9,19}, color_white},
    {{43,34,20}, {46,37,20}, {44,36,20}, color_white},
    {{45,33,21}, {18,7,21}, {46,37,21}, color_white},
    {{16,38,22}, {68,39,22}, {1,40,22}, color_white},
    {{42,31,23}, {17,4,23}, {45,33,23}, color_white},
    {{49,41,24}, {47,42,24}, {60,43,24}, color_white},
    {{54,44,25}, {53,45,25}, {48,46,25}, color_white},
    {{55,47,26}, {52,48,26}, {54,44,26}, color_white},
    {{58,49,27}, {56,50,27}, {50,51,27}, color_white},
    {{51,52,28}, {57,53,28}, {58,49,28}, color_white},
    {{47,42,29}, {64,54,29}, {60,43,29}, color_white},
    {{64,54,30}, {47,42,30}, {65,55,30}, color_white},
    {{49,56,31}, {59,57,31}, {67,58,31}, color_white},
    {{59,57,32}, {66,59,32}, {67,58,32}, color_white},
    {{78,60,33}, {77,61,33}, {5,2,33}, color_white},
    {{75,3,34}, {74,62,34}, {4,1,34}, color_white},
    {{73,63,35}, {72,64,35}, {3,65,35}, color_white},
    {{2,66,36}, {69,67,36}, {68,68,36}, color_white},
    {{83,69,37}, {82,70,37}, {84,71,37}, color_white},
    {{82,72,38}, {80,73,38}, {81,74,38}, color_white},
    {{85,75,39}, {84,71,39}, {86,76,39}, color_white},
    {{86,77,40}, {83,69,40}, {85,78,40}, color_white},
    {{80,79,41}, {83,80,41}, {86,76,41}, color_white},
    {{83,69,42}, {84,71,42}, {85,78,42}, color_white},
    {{84,71,43}, {81,81,43}, {86,76,43}, color_white},
    {{81,81,44}, {80,79,44}, {86,76,44}, color_white},
    {{82,82,45}, {81,81,45}, {84,71,45}, color_white},
    {{82,83,46}, {83,80,46}, {80,79,46}, color_white},
    {{87,84,47}, {91,85,47}, {92,86,47}, color_white},
    {{87,84,48}, {92,86,48}, {90,87,48}, color_white},
    {{92,88,49}, {89,89,49}, {90,90,49}, color_white},
    {{89,89,50}, {88,91,50}, {90,92,50}, color_white},
    {{88,91,51}, {87,93,51}, {90,94,51}, color_white},
    {{91,95,52}, {89,89,52}, {92,86,52}, color_white},
    {{91,96,53}, {88,91,53}, {89,89,53}, color_white},
    {{91,97,54}, {87,93,54}, {88,91,54}, color_white},
    {{6,98,55}, {70,99,55}, {7,100,55}, color_white},
    {{7,100,56}, {71,101,56}, {8,102,56}, color_white},
    {{9,103,57}, {71,101,57}, {72,64,57}, color_white},
    {{10,104,58}, {72,64,58}, {73,63,58}, color_white},
    {{10,104,59}, {74,62,59}, {11,105,59}, color_white},
    {{12,106,60}, {74,62,60}, {75,3,60}, color_white},
    {{12,106,61}, {76,107,61}, {13,108,61}, color_white},
    {{14,109,62}, {76,107,62}, {77,61,62}, color_white},
    {{15,110,63}, {77,61,63}, {78,60,63}, color_white},
    {{4,111,64}, {2,112,64}, {5,113,64}, color_white},
    {{15,110,65}, {79,114,65}, {16,38,65}, color_white},
    {{6,98,66}, {68,68,66}, {69,67,66}, color_white},
    {{11,115,67}, {15,116,67}, {7,117,67}, color_white},
    {{3,65,68}, {69,67,68}, {2,66,68}, color_white},
    {{32,18,69}, {27,118,69}, {28,22,69}, color_white},
    {{5,2,70}, {2,119,70}, {78,60,70}, color_white},
    {{3,65,71}, {4,1,71}, {73,63,71}, color_white},
    {{17,4,72}, {24,120,72}, {22,5,72}, color_white},
    {{18,7,73}, {19,6,73}, {23,8,73}, color_white},
    {{24,120,74}, {28,22,74}, {22,5,74}, color_white},
    {{19,6,75}, {22,5,75}, {23,8,75}, color_white},
    {{31,10,76}, {23,8,76}, {26,12,76}, color_white},
    {{22,5,77}, {28,22,77}, {25,11,77}, color_white},
    {{23,8,78}, {25,11,78}, {26,12,78}, color_white},
    {{24,120,79}, {21,20,79}, {28,22,79}, color_white},
    {{25,11,80}, {27,118,80}, {26,12,80}, color_white},
    {{25,11,81}, {28,22,81}, {27,118,81}, color_white},
    {{26,12,82}, {29,17,82}, {30,13,82}, color_white},
    {{33,19,83}, {30,13,83}, {29,17,83}, color_white},
    {{26,12,84}, {27,118,84}, {29,17,84}, color_white},
    {{28,22,85}, {35,21,85}, {32,18,85}, color_white},
    {{36,26,86}, {32,18,86}, {35,21,86}, color_white},
    {{33,19,87}, {36,26,87}, {34,14,87}, color_white},
    {{38,27,88}, {35,21,88}, {39,25,88}, color_white},
    {{34,14,89}, {36,26,89}, {41,28,89}, color_white},
    {{24,120,90}, {42,31,90}, {21,30,90}, color_white},
    {{34,14,91}, {41,28,91}, {37,15,91}, color_white},
    {{41,28,92}, {38,27,92}, {40,121,92}, color_white},
    {{38,27,93}, {39,25,93}, {40,121,93}, color_white},
    {{40,32,94}, {42,31,94}, {43,34,94}, color_white},
    {{41,122,95}, {43,34,95}, {44,36,95}, color_white},
    {{41,122,96}, {40,32,96}, {43,34,96}, color_white},
    {{44,36,97}, {46,37,97}, {20,9,97}, color_white},
    {{43,34,98}, {45,33,98}, {46,37,98}, color_white},
    {{45,33,99}, {17,4,99}, {18,7,99}, color_white},
    {{51,52,100}, {55,47,100}, {57,53,100}, color_white},
    {{46,37,101}, {18,7,101}, {20,9,101}, color_white},
    {{48,46,102}, {50,51,102}, {54,44,102}, color_white},
    {{50,51,103}, {48,123,103}, {58,49,103}, color_white},
    {{59,57,104}, {63,124,104}, {66,59,104}, color_white},
    {{48,46,105}, {53,45,105}, {51,125,105}, color_white},
    {{51,125,106}, {53,45,106}, {52,48,106}, color_white},
    {{53,45,107}, {54,44,107}, {52,48,107}, color_white},
    {{54,44,108}, {50,51,108}, {55,47,108}, color_white},
    {{50,51,109}, {56,50,109}, {55,47,109}, color_white},
    {{55,47,110}, {56,50,110}, {57,53,110}, color_white},
    {{56,50,111}, {58,49,111}, {57,53,111}, color_white},
    {{58,49,112}, {48,123,112}, {51,52,112}, color_white},
    {{65,55,113}, {47,42,113}, {49,56,113}, color_white},
    {{62,126,114}, {63,124,114}, {61,127,114}, color_white},
    {{49,41,115}, {60,43,115}, {59,128,115}, color_white},
    {{59,128,116}, {60,43,116}, {61,127,116}, color_white},
    {{60,43,117}, {62,126,117}, {61,127,117}, color_white},
    {{60,43,118}, {64,54,118}, {62,126,118}, color_white},
    {{63,124,119}, {65,55,119}, {66,59,119}, color_white},
    {{65,55,120}, {67,58,120}, {66,59,120}, color_white},
    {{65,55,121}, {49,56,121}, {67,58,121}, color_white},
    {{2,119,122}, {68,39,122}, {79,114,122}, color_white},
    {{2,119,123}, {79,114,123}, {78,60,123}, color_white},
    {{5,2,124}, {77,61,124}, {76,107,124}, color_white},
    {{5,2,125}, {76,107,125}, {75,3,125}, color_white},
    {{4,1,126}, {74,62,126}, {73,63,126}, color_white},
    {{3,65,127}, {72,64,127}, {71,101,127}, color_white},
    {{3,65,128}, {71,101,128}, {70,99,128}, color_white},
    {{16,38,129}, {79,114,129}, {68,39,129}, color_white},
    {{6,98,130}, {69,67,130}, {70,99,130}, color_white},
    {{7,100,131}, {70,99,131}, {71,101,131}, color_white},
    {{9,103,132}, {8,102,132}, {71,101,132}, color_white},
    {{10,104,133}, {9,103,133}, {72,64,133}, color_white},
    {{10,104,134}, {73,63,134}, {74,62,134}, color_white},
    {{12,106,135}, {11,105,135}, {74,62,135}, color_white},
    {{12,106,136}, {75,3,136}, {76,107,136}, color_white},
    {{14,109,137}, {13,108,137}, {76,107,137}, color_white},
    {{15,110,138}, {14,109,138}, {77,61,138}, color_white},
    {{4,111,139}, {3,129,139}, {2,112,139}, color_white},
    {{15,110,140}, {78,60,140}, {79,114,140}, color_white},
    {{6,98,141}, {1,130,141}, {68,68,141}, color_white},
    {{7,117,67}, {8,131,67}, {9,132,67}, color_white},
    {{9,132,67}, {10,133,67}, {7,117,67}, color_white},
    {{10,133,67}, {11,115,67}, {7,117,67}, color_white},
    {{11,115,67}, {12,134,67}, {13,135,67}, color_white},
    {{13,135,67}, {14,136,67}, {15,116,67}, color_white},
    {{15,116,67}, {16,137,67}, {1,138,67}, color_white},
    {{1,138,67}, {6,139,67}, {15,116,67}, color_white},
    {{6,139,67}, {7,117,67}, {15,116,67}, color_white},
    {{11,115,67}, {13,135,67}, {15,116,67}, color_white},
    {{3,65,142}, {70,99,142}, {69,67,142}, color_white},
    {{32,18,143}, {29,17,143}, {27,118,143}, color_white},
    {{33,19,144}, {34,14,144}, {30,13,144}, color_white},
    {{36,26,145}, {33,19,145}, {32,18,145}, color_white},
    {{38,27,146}, {36,26,146}, {35,21,146}, color_white},
    {{24,120,147}, {17,4,147}, {42,31,147}, color_white},
    {{51,52,148}, {52,48,148}, {55,47,148}, color_white},
    {{59,57,149}, {61,127,149}, {63,124,149}, color_white},
    {{62,126,150}, {64,54,150}, {63,124,150}, color_white},
    {{63,124,151}, {64,54,151}, {65,55,151}, color_white}}

function rabbitMesh()
   return meshdata(rabbitVertices, rabbitUVs, rabbitNormals, rabbitTriangles,
      -- backface
      true)
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

pico-8 cartridge // http://www.pico-8.com
version 29
__lua__
-- realtime renderering
-- by sthilaid

-->8
-- data and global vars
g_pal={0,128,133,5,134,6,135,15,7,3}
g_whiteShades={0,1,2,3,4,5,6,7,8} -- color indices for light shading from lighting 0->1
g_bgcolor=9
g_interpolateTriangles = true
g_performBackfaceCulling = true
g_renderMode = 2 -- 0: pixel draw 1: wire mesh 2: triangle draw

-- runtime globals
g_objects={}
g_lights={}
g_zbuff = {}
g_triangles = {}
g_camera=false
g_cube=false
g_apple=false
g_dirLight=false
g_pointLight=false
g_cubeRotMat=0
g_screenMat = 0
g_campitch=0
g_camyaw=0
g_camlen=3
g_visibleMeshIndex = 0

-- runtime debug
g_triCount = 0
g_totalTriCount = 0
g_pixelFragCount = 0
g_renderStartDT=0
g_renderEndDT=0
g_geometryDT=0
g_rasterDT=0
g_pixelDT=0

-------------------------------------------------------------------------------
-->8
-- flow

function _init()
   g_camera = makecam(g_campitch,g_camyaw,g_camlen)

   g_cube = object(mmult(transmatrix(vpoint(0,0,0)), rotmatrix(0.1,0.1,0)), cubemesh())
   g_cube.isVisible = false
   add(g_objects, g_cube)

   g_apple = object(rotmatrix(0.05,0.1,0), appleMesh())
   add(g_objects, g_apple)

   local plane = add(g_objects, object(matrix(), planemesh()))
   plane.isVisible =  false

   -- room = object(mmult(transmatrix(vpoint(0,0,0)), rotmatrix(0,0,0)), roommesh())
   -- add(g_objects, room)
   -- local rabbit = object(matrix(), rabbitMesh())
   -- add(g_objects, rabbit)

   g_dirLight = dirlight(vnormalize(vdir(-0.75,1,0)), 0.4)
   add(g_lights, g_dirLight)
   g_pointLight = pointlight(vpoint(-5,0,0), 4)
   add(g_lights, g_pointLight)

   g_cubeRotMat= rotmatrix(0.25/30, 0.1/30, 0.15/30)
   g_screenMat = mmult(transmatrix(vpoint(64, 64, 0.5)), scalematrix(64, 64, 0.5))
end

function _update()
   if (btn(⬆️)) g_camlen -= 0.1
   if (btn(⬇️)) g_camlen += 0.1
   if (btn(➡️)) g_camyaw += 0.01
   if (btn(⬅️)) g_camyaw -= 0.01
   if (btnp(❎)) g_renderMode = (g_renderMode+1) % 3
   if btnp(🅾️) then
      g_visibleMeshIndex = (g_visibleMeshIndex+1) % #g_objects + 1
      for i=1,#g_objects do
         local isVisible = false
         if (i == g_visibleMeshIndex) isVisible = true
         g_objects[i].isVisible = isVisible
      end
   end
   g_camera = makecam(g_campitch,g_camyaw,g_camlen)

   -- local translation = mgettranslation(cube.mat)
   -- mtranslate(g_cube.mat, vscale(translation, -1))
   -- g_cube.mat = mmult(g_cubeRotMat, cube.mat)
   -- mtranslate(cube.mat, translation)

   if (shouldRotateLight) g_dirLight.dir = mapply(rotmatrix(0,0.01,0), g_dirLight.dir)
   g_pointLight.pos = mapply(rotmatrix(0,0.01,0), g_pointLight.pos)
end

function _draw()
   for i,c in pairs(g_pal) do
      pal(i-1,c,1)
   end
   cls(g_bgcolor)
   render3d(g_camera, g_objects, g_lights)
   dflush()
end

-------------------------------------------------------------------------------
-- debug

debugStrs = {}
function dprint(str) add(debugStrs, str) end
function dflush()
   color(15)
   print("mem:"..(stat(0)/2048).."%".." cpu:"..((g_renderEndDT-g_renderStartDT) * 100).."% @"..stat(7).."fps")
   print("tri: "..g_triCount.."/"..g_totalTriCount.." pfrags: "..g_pixelFragCount)
   print("g:"..g_geometryDT.." r:"..g_rasterDT.." p:"..g_pixelDT)
   for s in all(debugStrs) do
      print(s)
   end

   local modeStr = "triangle draw"
   if (g_renderMode == 0) modeStr = "pixel draw"
   if (g_renderMode == 1) modeStr = "wiremesh draw"
   print("mode: "..modeStr, 0, 116)
   print("⬅️➡️⬆️⬇️:cam ❎:mode 🅾️:mesh", 0, 122)
   debugStrs = {}
   g_geometryDT, g_rasterDT, g_pixelDT, g_triCount, g_totalTriCount, g_pixelFragCount = 0,0,0,0,0,0
end

-->8
-- math

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
function abs(x) if x<0 then return -x else return x end end
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
function vsub(v1, v2) return {v1[1]-v2[1], v1[2]-v2[2], v1[3]-v2[3], v1[4]+v2[4]} end
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

function matrix() return {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} end
function rotmatrix(x,y,z) -- ammount to rotate around each axis
   function rot_x(xangle)
      local m = matrix()
      m[6] = cos(xangle)
      m[7] = -sin(xangle)
      m[10] = sin(xangle)
      m[11] = cos(xangle)
      return m
   end
   function rot_y(yangle)
      local m = matrix()
      m[1] = cos(yangle)
      m[3] = sin(yangle)
      m[9] = -sin(yangle)
      m[11] = cos(yangle)
      return m
   end
   function rot_z(zangle)
      local m = matrix()
      m[1] = cos(zangle)
      m[2] = -sin(zangle)
      m[5] = sin(zangle)
      m[6] = cos(yangle)
      return m
   end

   return mmult(mmult(rot_x(x), rot_y(y)), rot_z(z))
end
function transmatrix(t)
   local m = matrix()
   m[4] = t[1]
   m[8] = t[2]
   m[12]= t[3]
   return m
end
function scalematrix(sx, sy, sz)
   local m = matrix()
   m[1] = sx
   m[6] = sy
   m[11] = sz
   m[16] = 1
   return m
end
function mgettranslation(m)
   return vpoint(m[4], m[8], m[12])
end
function mtranslate(m, t)
   m[4] += t[1]
   m[8] += t[2]
   m[12] += t[3]
end
function ortho(l, r, t, b, n, f)
   local m = matrix()
   m[1] = 2 / (r - l)
   m[6] = 2 / (t - b)
   m[11] = 2 / (f - n)
   m[4] = -(r+l) / (r-l)
   m[8] = -(t+b) / (t-b)
   m[12] = -(f+n) / (f-n)
   m[16] = 1
   return m
end
function perspective(l, r, t, b, n, f)
   local m = matrix()
   m[1] = 2*n / (r - l)
   m[6] = 2*n / (t - b)
   m[11] = (f+n) / (f-n)
   m[3] = -(r+l) / (r-l)
   m[7] = -(t+b) / (t-b)
   m[12] = -2*f*n / (f-n)
   m[15] = 1
   m[16] = 0
   return m
end
-- function perspective2(fovy, aspect, near, far)
--    local theta = fovy/2
--    local top = near * sin(theta) / cos(theta) -- tan(x) = sin(x) / cos(x)
--    local bottom = -top
--    local right = top * aspect
--    local left = -right
--    return perspective(left, right, bottom, top, near, far);
-- end
function mmult(m1, m2)
   local m = matrix()
   for i=1,4 do
      for j=1,4 do
         local m_ij = 0
         for r=1,4 do
            m_ij += m1[(i-1)*4+r] * m2[(r-1)*4+j]
         end
         m[(i-1)*4+j] = m_ij
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
         res[i] += m[(i-1)*4+j] * v[j]
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

function meshdata(vertices, normals, triangles, cullBackface)
   return {verts=vertices,
           normals=normals,
           tris=triangles,
           cullBackface=cullBackface}
end
function object(mat, mesh) return {mat=mat, mesh=mesh, isVisible=true} end

function makecam(pitch,yaw, len)
   return mmult(mmult(perspective(-3, 3, 3, -3, -2, -100), transmatrix(vpoint(0,0,len))), rotmatrix(pitch,yaw,0))
end

-- lights: 0-directional, 1-point
function dirlight(dir, r0) return {type=0, pos=vpoint(0,0,0), dir=dir, r0=r0} end
function pointlight(pos, r0) return {type=1, pos=pos, dir=vdir(0,0,0), r0=r0} end
function lightdata(ltype,ldir,lr0,lr) return {type=ltype, dir=ldir, r0=lr0, r=lr} end

function fragment(v, n)
   return {vx=v[1], vy=v[2], vz=v[3],
           nx=n[1], ny=n[2], nz=n[3],
           lightsData = false,
   }
end

-------------------------------------------------------------------------------
-->8
-- rendering pipeline

function geometryVertexShading(projectionMatrix, mesh, v_index, n_index,
                               processedVerticesCache, lights, tcolor)
   local v = mesh.verts[v_index]

   -- lighting
   local lightsData  = {}

   for l in all(lights) do
      local ltype,ldir,lr0,lr = false,false,false,false
      if l then
         if l.type == 0 then ldir = l.dir
         elseif l.type == 1 then
            local deltaToLight = vsub(l.pos, v)
            lr    = max(0.001, vnorm(deltaToLight))
            ldir  = vscale(deltaToLight, 1/lr)
         end
         ltype= l.type
         lr0  = l.r0
         add(lightsData, lightdata(ltype,ldir,lr0,lr))
      end
   end
   
   -- vertex projection
   local camspaceVert = processedVerticesCache[v_index]
   if not camspaceVert then
      camspaceVert  = mapply(projectionMatrix, v)
      camspaceVert  = vscale(camspaceVert, 1/camspaceVert[4]) // normalize the w component
   end
   local n  = mesh.normals[n_index]

   -- fragment creation
   local vfrag  = fragment(camspaceVert, n)
   vfrag.lightsData  = lightsData
   return vfrag
end

function geometryClipping(projVertFrag, clippedTriangles)
    -- todo
   add(clippedTriangles, projVertFrag)
end

function processGeometries(cam, objects, lights)
   for obj in all(objects) do
      if (not obj.isVisible) goto nextobj
      local mesh = obj.mesh
      local processedVerticesCache = {}
      local projectionMatrix = mmult(g_screenMat, mmult(cam, obj.mat))
      for t in all(mesh.tris) do
         local startDT = stat(1)
         local vfrags= {}
         local clippedTriangleFrags = {}
         local triangleColor = t[4]
         for tindex=1,3 do
            local v_index   = t[tindex][1]
            local n_index   = t[tindex][2]
            local vFrag     = geometryVertexShading(projectionMatrix, mesh, v_index, n_index,
                                                    processedVerticesCache, lights, triangleColor)
            --dprint("pv: "..vstr(vFrag.n))
            add(vfrags, vFrag)
         end
         
         -- backface culling
         if (not g_performBackfaceCulling) or triangleSignedArea(vfrags) < 0 then
            g_triCount += 1
            g_geometryDT += stat(1) - startDT

            if (g_renderMode == 0) rasterizeTriangle(vfrags)
            if (g_renderMode == 1) rasterizeTriangleToEdges(vfrags)
            if (g_renderMode == 2) storeTriangle(vfrags)
         else
            g_geometryDT += stat(1) - startDT
         end
         g_totalTriCount += 1
      end
      ::nextobj::
   end
   if g_renderMode == 2 then -- triangle draw mode
      qsort(g_triangles,1,#g_triangles,'depth',lt)
      for t in all(g_triangles) do
         drawTriangle(t)
      end
   end
end

function triangleSignedArea(vfrags)
   -- sign of the z component of the cross product of the edge vectors (p2-p1 and p3-p1)
   local p1x, p1y = vfrags[1].vx, vfrags[1].vy
   local p2x, p2y = vfrags[2].vx, vfrags[2].vy
   local p3x, p3y = vfrags[3].vx, vfrags[3].vy
   return -p1y*p2x+p1x*p2y+p1y*p3x-p2y*p3x-p1x*p3y+p2x*p3y
end
function computeEdgeParams(p1x,p1y,p2x,p2y) return {a=-(p2y-p1y), b=(p2x-p1x), c=(p2y-p1y)*p1x - (p2x-p1x)*p1y} end
function edgeSign(px,py,a,b,c) return a*px + b*py + c end
function baryInterpVertex(w,u,v,v1,v2,v3) return vadd(vadd(vscale(v1, w), vscale(v2, u)), vscale(v3, v)) end
function baryInterpPoint(w,u,v,v1,v2,v3) return padd(padd(pscale(v1, w), pscale(v2, u)), pscale(v3, v)) end
function baryInterpValue(w,u,v,v1,v2,v3) return v1*w + v2*u + v3*v end
function interpolateTri(x, y, vfrag1, vfrag2, vfrag3, edgeValues)
   local areaSum = (edgeValues[1]  + edgeValues[2] + edgeValues[3])
   -- should use perspective corrected coordinates?
   local barycentric_u = edgeValues[2] / areaSum
   local barycentric_v = edgeValues[3] / areaSum
   local barycentric_w = 1 - barycentric_u - barycentric_v
   -- only interpolate depth instead of whole vertex pos?
   local depth  = baryInterpValue(barycentric_w, barycentric_u, barycentric_v, vfrag1.vz, vfrag2.vz, vfrag3.vz)
   local n      = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v,
                                   vdir(vfrag1.nx, vfrag1.ny, vfrag1.nz),
                                   vdir(vfrag2.nx, vfrag2.ny, vfrag2.nz),
                                   vdir(vfrag3.nx, vfrag3.ny, vfrag3.nz))

   local interpLightsData = {}
   for lindex=1,#vfrag1.lightsData do
      vldata1 = vfrag1.lightsData[lindex]
      vldata2 = vfrag2.lightsData[lindex]
      vldata3 = vfrag3.lightsData[lindex]
      local ltype,ldir,lr0,lr = vldata1.type,false,false,false
      if ltype == 0 then -- directional
         ldir   = vldata1.dir
         lr0    = vldata1.r0
      elseif ltype == 1 then -- omni
         -- ldir   = baryInterpVertex(barycentric_w, barycentric_u, barycentric_v,
         --                           vldata1.dir, vldata2.dir, vldata3.dir)
         ldir   = vldata1.dir
         lr0    = vldata1.r0
         lr     = baryInterpValue(barycentric_w, barycentric_u, barycentric_v,
                                  vldata1.r, vldata2.r, vldata3.r)
      end
      add(interpLightsData, lightdata(ltype,ldir,lr0,lr))
   end
   local pfrag  = fragment(vpoint(x,y,depth), n)
   pfrag.lightsData  = interpLightsData
   return pfrag
end

function rasterizeTriangle(vfrags)
   local p1x,p1y = vfrags[1].vx, vfrags[1].vy
   local p2x,p2y = vfrags[2].vx, vfrags[2].vy
   local p3x,p3y = vfrags[3].vx, vfrags[3].vy
   local edgeParams = {computeEdgeParams(p1x,p1y, p2x,p2y),
                       computeEdgeParams(p2x,p2y, p3x,p3y),
                       computeEdgeParams(p3x,p3y, p1x,p1y)}
   local pmin_x = flr(min(min(p1x, p2x), p3x))
   local pmin_y = flr(min(min(p1y, p2y), p3y))
   local pmax_x = flr(max(max(p1x, p2x), p3x))
   local pmax_y = flr(max(max(p1y, p2y), p3y))
   --dprint("min: "..pmin_x..","..pmin_y.." max: "..pmax_x..","..pmax_y)
   --dprint("params: a:"..edgeParams[1].a.."b: "..edgeParams[1].b.."c: "..edgeParams[1].c)
   for x=pmin_x, pmax_x do
      for y=pmin_y, pmax_y do
         local pixelStartDT = stat(1)
         local isInside = true
         local edgeValues = {}
         for e=1,3 do
            edgeValues[e] = edgeSign(x,y, edgeParams[e].a, edgeParams[e].b, edgeParams[e].c)
            -- can use this instead if we don't want backface culling at the pixel level
            --isInside = e==1 or edgeValues[e] == 0 or sign(edgeValues[e]) == sign(edgeValues[e-1])
            isInside = isInside and edgeValues[e] <= 0
            if not isInside then
               break
            end
         end
         --dprint(edgeValues[1])
         --dprint("("..x..","..y.."): "..edgeValues[1]..","..edgeValues[2]..","..edgeValues[3])
         if isInside then
            local pfrag = false
            if g_interpolateTriangles then
               pfrag = interpolateTri(x, y, vfrags[1], vfrags[2], vfrags[3], edgeValues)
            else
               local v = vpoint(x, y, vfrags[1].vz)
               local n = vdir(vfrags[1].nx, vfrags[1].ny, vfrags[1].nz)
               pfrag = fragment(v, n)
               pfrag.lightsData  = vfrags[1].lightsData
            end

            g_rasterDT += stat(1) - pixelStartDT
            processPixelFragment(pfrag)
         end
      end
   end
end

function rasterizeEdge(p1x,p1y,p1z,p2x,p2y,p2z)
   line(p1x,p1y,p2x,p2y,8)
end
function rasterizeTriangleToEdges(vfrags)
   local p1x,p1y,p1z = vfrags[1].vx, vfrags[1].vy, vfrags[1].vz
   local p2x,p2y,p2z = vfrags[2].vx, vfrags[2].vy, vfrags[2].vz
   local p3x,p3y,p3z = vfrags[3].vx, vfrags[3].vy, vfrags[3].vz
   rasterizeEdge(p1x,p1y,p1z,p2x,p2y,p2z)
   rasterizeEdge(p2x,p2y,p2z,p3x,p3y,p3z)
   rasterizeEdge(p3x,p3y,p3z,p1x,p1y,p1z)
end

function getShadedColor(lightsData, nx, ny, nz)
   local lightRatio = 0.0
   for ldata in all(lightsData) do
      if ldata.type == 0 then
         lightRatio += clamp(0,1,max(0,ldata.r0 * vdot(ldata.dir, vdir(nx,ny,nz))))
      elseif ldata.type == 1 then
         lightRatio += clamp(0,1,max(0, vdot(ldata.dir, vdir(nx,ny,nz)) * ldata.r0*ldata.r0 / max(ldata.r*ldata.r, 0.0001)))
      end
   end
   lightRatio = clamp(0,1,lightRatio)
   local shadedColorIndex = min(flr(lightRatio * (#g_whiteShades-1))+1, #g_whiteShades)
   local shadedColor = g_whiteShades[shadedColorIndex]
   return shadedColor
end

function buffGetValue(buf, x, y) return buf[y*128+x] end
function buffSetValue(buf, x, y, val) buf[y*128+x] = val end
function processPixelFragment(pfrag)
   local startDT = stat(1)
   g_pixelFragCount += 1

   local currentZValue = buffGetValue(g_zbuff, pfrag.vx, pfrag.vy)
   if (currentZValue and pfrag.vz < currentZValue) return
   buffSetValue(g_zbuff, pfrag.vx, pfrag.vy, pfrag.vz)
   
   -- lighting
   local shadedColor = getShadedColor(pfrag.lightsData, pfrag.nx,pfrag.ny, pfrag.nz)
   --print(shadedColor)
   
   -- pixel render
   pset(pfrag.vx, pfrag.vy, shadedColor)
   g_pixelDT += stat(1) - startDT
end

function storeTriangle(vfrags)
   vfrags.depth = max(max(vfrags[1].vz, vfrags[2].vz), vfrags[3].vz)
   add(g_triangles,vfrags)
end

function drawTriangle(vfrags)
   local p1x,p1y,p1z = vfrags[1].vx, vfrags[1].vy, vfrags[1].vz
   local p2x,p2y,p2z = vfrags[2].vx, vfrags[2].vy, vfrags[2].vz
   local p3x,p3y,p3z = vfrags[3].vx, vfrags[3].vy, vfrags[3].vz

   local pfrag=vfrags[1]
   local shadedColor = getShadedColor(pfrag.lightsData, pfrag.nx,pfrag.ny,pfrag.nz)
   trifill(p1x,p1y,p2x,p2y,p3x,p3y,shadedColor)
end

function render3d(camera, objects, lights)
   g_renderStartDT = stat(1)
   g_zbuff, g_triangles = {}, {}
   processGeometries(camera, objects, lights)

   g_renderEndDT = stat(1)
end

-------------------------------------------------------------------------------
-->8
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
      -- normals
      {vdir(0.0000, 1.0000, 0.0000),
       vdir(0.0000, 0.0000, 1.0000),
       vdir(-1.0000, 0.0000, 0.0000),
       vdir(0.0000, -1.0000, 0.0000),
       vdir(1.0000, 0.0000, 0.0000),
       vdir(0.0000, 0.0000, -1.0000)},
      -- triangles (v,n)
      {  {{5,1}, {3,1}, {1,1}},
         {{3,2}, {8,2}, {4,2}},
         {{7,3}, {6,3}, {8,3}},
         {{2,4}, {8,4}, {6,4}},
         {{1,5}, {4,5}, {2,5}},
         {{5,6}, {2,6}, {6,6}},
         {{5,1}, {7,1}, {3,1}},
         {{3,2}, {7,2}, {8,2}},
         {{7,3}, {5,3}, {6,3}},
         {{2,4}, {4,4}, {8,4}},
         {{1,5}, {3,5}, {4,5}},
         {{5,6}, {1,6}, {2,6}}})
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
        {vdir(0.0000, 0.0000, -1.0000),
         vdir(1.0000, 0.0000, 0.0000),
         vdir(0.0000, 1.0000, 0.0000),
         vdir(-1.0000, 0.0000, 0.0000)},
        {{{8,1}, {3,1}, {4,1}},
         {{6,2}, {7,2}, {8,2}},
         {{8,3}, {2,3}, {6,3}},
         {{4,4}, {1,4}, {2,4}},
         {{8,1}, {7,1}, {3,1}},
         {{6,2}, {5,2}, {7,2}},
         {{8,3}, {4,3}, {2,3}},
         {{4,4}, {3,4}, {1,4}}})
end
function planemesh()
   return meshdata(
        {vpoint(1.000000, 1.000000, 0.000000),
         vpoint(1.000000, -1.000000, 0.000000),
         vpoint(-1.000000, -1.000000, 0.000000),
         vpoint(-1.000000, 1.000000, 0.000000)},
        {vdir(0.0000, 0.0000, -1.0000)},
        {{{1,1}, {2,1}, {3,1}},
         {{1,1}, {3,1}, {4,1}}})
end
appleVertices = {
    vpoint(0.005404, 0.761771, -0.500399),
    vpoint(-0.010084, 0.482872, -0.839171),
    vpoint(-0.012526, -0.000000, -0.968364),
    vpoint(0.023313, -0.458409, -0.844697),
    vpoint(0.035651, -0.692655, -0.470820),
    vpoint(0.005190, 0.657124, -0.000584),
    vpoint(0.431062, 0.782246, -0.253032),
    vpoint(0.686725, 0.437787, -0.401163),
    vpoint(0.864982, 0.006544, -0.501204),
    vpoint(0.771063, -0.446717, -0.431835),
    vpoint(0.433013, -0.717226, -0.250000),
    vpoint(0.423549, 0.753163, 0.247459),
    vpoint(0.692789, 0.455101, 0.426907),
    vpoint(0.861803, -0.003314, 0.499549),
    vpoint(0.750000, -0.500000, 0.433013),
    vpoint(0.460420, -0.674564, 0.241328),
    vpoint(-0.000000, 0.756350, 0.500000),
    vpoint(-0.000000, 0.500000, 0.866025),
    vpoint(-0.000000, -0.000000, 1.000000),
    vpoint(0.002474, -0.425387, 0.889564),
    vpoint(-0.008560, -0.629215, 0.528972),
    vpoint(0.000000, -0.594428, -0.000000),
    vpoint(-0.429401, 0.756350, 0.250900),
    vpoint(-0.746800, 0.501590, 0.433521),
    vpoint(-0.802415, 0.073784, 0.489289),
    vpoint(-0.746800, -0.498410, 0.433521),
    vpoint(-0.433013, -0.717226, 0.250000),
    vpoint(-0.385427, 0.751821, -0.237312),
    vpoint(-0.672266, 0.476368, -0.409337),
    vpoint(-0.860945, -0.001744, -0.498417),
    vpoint(-0.750000, -0.500000, -0.433013),
    vpoint(-0.433013, -0.717226, -0.250000),
    vpoint(0.005190, 0.959313, -0.050042),
    vpoint(0.048022, 0.959313, -0.025313),
    vpoint(0.005190, 0.672644, -0.050042),
    vpoint(0.005190, 0.959313, -0.050042),
    vpoint(0.048022, 0.672644, -0.025313),
    vpoint(0.048022, 0.959313, -0.025313),
    vpoint(0.048022, 0.672644, 0.024145),
    vpoint(0.048022, 0.959313, 0.024145),
    vpoint(0.005190, 0.672644, 0.048874),
    vpoint(0.005190, 0.959313, 0.048874),
    vpoint(-0.037642, 0.672644, 0.024145),
    vpoint(-0.037642, 0.959313, 0.024145),
    vpoint(-0.037642, 0.672644, -0.025313),
    vpoint(-0.037642, 0.959313, -0.025313),
    vpoint(0.048022, 0.959313, 0.024145),
    vpoint(0.005190, 0.959313, 0.048874),
    vpoint(-0.037642, 0.959313, 0.024145),
    vpoint(-0.037642, 0.959313, -0.025313),
    vpoint(-0.034500, 1.184826, -0.028135),
    vpoint(-0.076851, 1.171124, -0.006577),
    vpoint(-0.078685, 1.164615, 0.042416),
    vpoint(-0.038168, 1.171809, 0.069852),
    vpoint(0.004183, 1.185512, 0.048295),
    vpoint(0.006017, 1.192021, -0.000699)
}

appleNormals = {
    vdir(-0.1651, -0.9679, 0.1894),
    vdir(0.4733, -0.2283, -0.8508),
    vdir(0.4153, 0.6929, -0.5894),
    vdir(0.2833, -0.7532, -0.5937),
    vdir(0.4344, 0.3701, -0.8212),
    vdir(-0.1639, 0.9656, 0.2021),
    vdir(0.9752, -0.2209, 0.0101),
    vdir(0.7477, 0.6637, -0.0194),
    vdir(0.6348, -0.7720, 0.0316),
    vdir(0.9247, 0.3807, 0.0067),
    vdir(-0.2528, 0.9661, 0.0523),
    vdir(-0.2196, -0.9708, 0.0965),
    vdir(0.3846, 0.7561, 0.5295),
    vdir(0.2057, -0.8580, 0.4706),
    vdir(0.4788, 0.3069, 0.8225),
    vdir(-0.1086, 0.9749, -0.1944),
    vdir(-0.1368, -0.9883, -0.0672),
    vdir(0.4917, -0.2234, 0.8416),
    vdir(-0.2143, -0.7896, 0.5751),
    vdir(-0.5088, 0.2228, 0.8316),
    vdir(0.1114, 0.9750, -0.1921),
    vdir(0.2402, -0.9689, -0.0598),
    vdir(-0.5059, -0.1322, 0.8524),
    vdir(-0.3138, 0.7777, 0.5447),
    vdir(-0.5653, -0.8249, -0.0000),
    vdir(-0.9316, 0.3625, 0.0275),
    vdir(0.2289, 0.9734, 0.0116),
    vdir(0.2728, -0.9621, -0.0000),
    vdir(-0.9762, -0.2168, 0.0040),
    vdir(-0.6516, 0.7543, -0.0802),
    vdir(-0.1992, -0.8334, -0.5156),
    vdir(-0.4572, 0.3340, -0.8242),
    vdir(0.1123, 0.9726, 0.2037),
    vdir(0.1509, -0.9654, 0.2128),
    vdir(-0.4419, -0.2657, -0.8568),
    vdir(-0.3816, 0.7221, -0.5770),
    vdir(0.0000, 0.0000, 1.0000),
    vdir(0.5000, 0.0000, -0.8660),
    vdir(1.0000, 0.0000, 0.0000),
    vdir(0.5000, 0.0000, 0.8660),
    vdir(-0.5000, 0.0000, 0.8660),
    vdir(0.9811, 0.1836, 0.0611),
    vdir(-1.0000, 0.0000, 0.0000),
    vdir(-0.5000, 0.0000, -0.8660),
    vdir(0.0000, -1.0000, 0.0000),
    vdir(-0.2507, 0.9608, 0.1182),
    vdir(0.4919, 0.1789, -0.8521),
    vdir(-0.5000, -0.0039, -0.8660),
    vdir(-0.5295, -0.1796, 0.8291),
    vdir(-0.9833, -0.1820, 0.0000),
    vdir(0.4525, 0.0043, 0.8917),
    vdir(0.4618, -0.1972, -0.8648),
    vdir(0.3841, 0.5904, -0.7099),
    vdir(0.2964, -0.8049, -0.5141),
    vdir(0.5295, 0.2167, -0.8201),
    vdir(0.9792, -0.2027, 0.0011),
    vdir(0.8115, 0.5825, 0.0460),
    vdir(0.5347, -0.8442, -0.0390),
    vdir(0.9389, 0.3439, -0.0141),
    vdir(0.3692, 0.7009, 0.6103),
    vdir(0.2059, -0.8546, 0.4768),
    vdir(0.5324, 0.2191, 0.8176),
    vdir(0.4896, -0.2164, 0.8446),
    vdir(-0.2130, -0.8478, 0.4857),
    vdir(-0.4931, 0.1752, 0.8521),
    vdir(-0.5384, -0.2147, 0.8149),
    vdir(-0.3142, 0.7797, 0.5416),
    vdir(-0.5706, -0.8212, 0.0036),
    vdir(-0.9890, 0.1167, -0.0909),
    vdir(-0.9924, -0.1030, 0.0667),
    vdir(-0.6698, 0.7395, -0.0672),
    vdir(-0.2428, -0.8074, -0.5377),
    vdir(-0.5325, 0.2213, -0.8170),
    vdir(-0.4726, -0.2173, -0.8540),
    vdir(-0.3864, 0.7434, -0.5460),
    vdir(0.9841, 0.1776, 0.0000),
    vdir(-0.2507, 0.9608, 0.1183),
    vdir(0.5308, 0.1740, -0.8295),
    vdir(-0.4524, -0.0049, -0.8918),
    vdir(-0.4914, -0.1843, 0.8512),
    vdir(-0.9798, -0.1904, -0.0620),
    vdir(0.5000, 0.0044, 0.8660)
}

appleTriangles = {
{{22,1}, {5,1}, {11,1}},
{{4,2}, {9,2}, {10,2}},
{{1,3}, {8,3}, {2,3}},
{{5,4}, {10,4}, {11,4}},
{{3,5}, {8,5}, {9,5}},
{{1,6}, {6,6}, {7,6}},
{{10,7}, {14,7}, {15,7}},
{{8,8}, {12,8}, {13,8}},
{{10,9}, {16,9}, {11,9}},
{{8,10}, {14,10}, {9,10}},
{{7,11}, {6,11}, {12,11}},
{{22,12}, {11,12}, {16,12}},
{{13,13}, {17,13}, {18,13}},
{{15,14}, {21,14}, {16,14}},
{{13,15}, {19,15}, {14,15}},
{{12,16}, {6,16}, {17,16}},
{{22,17}, {16,17}, {21,17}},
{{14,18}, {20,18}, {15,18}},
{{21,19}, {26,19}, {27,19}},
{{18,20}, {25,20}, {19,20}},
{{17,21}, {6,21}, {23,21}},
{{22,22}, {21,22}, {27,22}},
{{20,23}, {25,23}, {26,23}},
{{17,24}, {24,24}, {18,24}},
{{27,25}, {31,25}, {32,25}},
{{25,26}, {29,26}, {30,26}},
{{23,27}, {6,27}, {28,27}},
{{22,28}, {27,28}, {32,28}},
{{26,29}, {30,29}, {31,29}},
{{23,30}, {29,30}, {24,30}},
{{32,31}, {4,31}, {5,31}},
{{29,32}, {3,32}, {30,32}},
{{28,33}, {6,33}, {1,33}},
{{22,34}, {32,34}, {5,34}},
{{31,35}, {3,35}, {4,35}},
{{29,36}, {1,36}, {2,36}},
{{46,37}, {49,37}, {50,37}},
{{36,38}, {37,38}, {35,38}},
{{44,37}, {48,37}, {49,37}},
{{38,39}, {39,39}, {37,39}},
{{42,37}, {47,37}, {48,37}},
{{40,40}, {41,40}, {39,40}},
{{38,37}, {33,37}, {34,37}},
{{42,41}, {43,41}, {41,41}},
{{36,37}, {50,37}, {33,37}},
{{47,42}, {56,42}, {55,42}},
{{44,43}, {45,43}, {43,43}},
{{40,37}, {34,37}, {47,37}},
{{46,44}, {35,44}, {45,44}},
{{41,45}, {43,45}, {45,45}},
{{53,46}, {54,46}, {55,46}},
{{33,47}, {56,47}, {34,47}},
{{50,48}, {51,48}, {33,48}},
{{49,49}, {54,49}, {53,49}},
{{49,50}, {52,50}, {50,50}},
{{48,51}, {55,51}, {54,51}},
{{4,52}, {3,52}, {9,52}},
{{1,53}, {7,53}, {8,53}},
{{5,54}, {4,54}, {10,54}},
{{3,55}, {2,55}, {8,55}},
{{10,56}, {9,56}, {14,56}},
{{8,57}, {7,57}, {12,57}},
{{10,58}, {15,58}, {16,58}},
{{8,59}, {13,59}, {14,59}},
{{13,60}, {12,60}, {17,60}},
{{15,61}, {20,61}, {21,61}},
{{13,62}, {18,62}, {19,62}},
{{14,63}, {19,63}, {20,63}},
{{21,64}, {20,64}, {26,64}},
{{18,65}, {24,65}, {25,65}},
{{20,66}, {19,66}, {25,66}},
{{17,67}, {23,67}, {24,67}},
{{27,68}, {26,68}, {31,68}},
{{25,69}, {24,69}, {29,69}},
{{26,70}, {25,70}, {30,70}},
{{23,71}, {28,71}, {29,71}},
{{32,72}, {31,72}, {4,72}},
{{29,73}, {2,73}, {3,73}},
{{31,74}, {30,74}, {3,74}},
{{29,75}, {28,75}, {1,75}},
{{46,37}, {44,37}, {49,37}},
{{36,38}, {38,38}, {37,38}},
{{44,37}, {42,37}, {48,37}},
{{38,39}, {40,39}, {39,39}},
{{42,37}, {40,37}, {47,37}},
{{40,40}, {42,40}, {41,40}},
{{38,37}, {36,37}, {33,37}},
{{42,41}, {44,41}, {43,41}},
{{36,37}, {46,37}, {50,37}},
{{47,76}, {34,76}, {56,76}},
{{44,43}, {46,43}, {45,43}},
{{40,37}, {38,37}, {34,37}},
{{46,44}, {36,44}, {35,44}},
{{45,45}, {35,45}, {41,45}},
{{35,45}, {37,45}, {41,45}},
{{37,45}, {39,45}, {41,45}},
{{55,77}, {56,77}, {51,77}},
{{51,77}, {52,77}, {55,77}},
{{52,77}, {53,77}, {55,77}},
{{33,78}, {51,78}, {56,78}},
{{50,79}, {52,79}, {51,79}},
{{49,80}, {48,80}, {54,80}},
{{49,81}, {53,81}, {52,81}},
{{48,82}, {47,82}, {55,82}},
}

function appleMesh()
   return meshdata(appleVertices, appleNormals, appleTriangles)
end

-- rabbitVertices = {
--    vpoint(0.000000, -1.391679, -0.797764),
--    vpoint(-0.135392, 0.611448, -0.507234),
--    vpoint(0.517557, 0.609942, 0.052065),
--    vpoint(0.056643, 0.604025, 0.500860),
--    vpoint(-0.473639, 0.610471, 0.213351),
--    vpoint(0.398882, -1.391679, -0.690884),
--    vpoint(0.690884, -1.391679, -0.398882),
--    vpoint(0.797764, -1.391679, 0.000000),
--    vpoint(0.690884, -1.391679, 0.398882),
--    vpoint(0.398882, -1.391679, 0.690884),
--    vpoint(0.000000, -1.391679, 0.797764),
--    vpoint(-0.398882, -1.391679, 0.690884),
--    vpoint(-0.690883, -1.391679, 0.398882),
--    vpoint(-0.797764, -1.391679, 0.000000),
--    vpoint(-0.690884, -1.391679, -0.398881),
--    vpoint(-0.398882, -1.391679, -0.690883),
--    vpoint(0.007565, 1.280805, -0.534633),
--    vpoint(-0.001319, 0.824538, -0.592229),
--    vpoint(0.310472, 0.972387, -0.526726),
--    vpoint(-0.004148, 0.435775, -0.283392),
--    vpoint(-0.085475, 1.602911, -0.022652),
--    vpoint(0.457809, 1.236295, -0.356192),
--    vpoint(0.386913, 0.649667, -0.355688),
--    vpoint(0.149670, 1.545669, -0.219466),
--    vpoint(0.619546, 1.066289, 0.001204),
--    vpoint(0.504853, 0.650576, 0.163492),
--    vpoint(0.517314, 1.077296, 0.298671),
--    vpoint(0.330113, 1.505936, 0.119562),
--    vpoint(0.300922, 0.995839, 0.532239),
--    vpoint(0.246869, 0.529812, 0.342540),
--    vpoint(0.157971, 0.390043, -0.016812),
--    vpoint(0.129921, 1.364891, 0.469963),
--    vpoint(-0.000000, 0.963595, 0.612501),
--    vpoint(-0.133589, 0.587250, 0.454342),
--    vpoint(-0.162013, 1.455557, 0.363314),
--    vpoint(-0.306251, 1.004631, 0.530442),
--    vpoint(-0.117633, 0.423134, 0.210932),
--    vpoint(-0.519327, 1.063339, 0.300116),
--    vpoint(-0.402400, 1.433085, 0.146065),
--    vpoint(-0.615259, 1.078989, 0.007206),
--    vpoint(-0.497713, 0.646425, 0.157127),
--    vpoint(-0.326746, 1.373354, -0.349210),
--    vpoint(-0.537144, 0.901938, -0.308731),
--    vpoint(-0.267728, 0.433814, -0.075064),
--    vpoint(-0.299692, 1.063339, -0.519572),
--    vpoint(-0.273868, 0.687422, -0.459224),
--    vpoint(0.403289, 2.117685, 0.035751),
--    vpoint(-0.242199, 2.117975, -0.073327),
--    vpoint(0.212305, 2.118290, -0.061642),
--    vpoint(-0.357801, 2.117975, 0.073328),
--    vpoint(-0.368098, 1.280821, -0.110259),
--    vpoint(-0.102704, 1.264386, -0.020220),
--    vpoint(-0.095093, 1.641306, -0.080582),
--    vpoint(-0.045583, 1.810991, 0.025934),
--    vpoint(-0.231820, 1.283460, 0.110408),
--    vpoint(-0.504907, 1.641307, 0.080582),
--    vpoint(-0.497220, 1.264120, 0.020212),
--    vpoint(-0.554417, 1.810991, -0.025934),
--    vpoint(0.231902, 1.280821, -0.110259),
--    vpoint(0.532326, 1.829810, -0.061515),
--    vpoint(0.497296, 1.264386, -0.020220),
--    vpoint(0.552316, 1.641104, 0.024805),
--    vpoint(0.354862, 1.266797, 0.089824),
--    vpoint(0.422333, 1.878187, 0.116543),
--    vpoint(0.073201, 1.844584, 0.066144),
--    vpoint(0.102780, 1.264120, 0.020212),
--    vpoint(0.047684, 1.641104, -0.024805),
--    vpoint(0.007158, 0.302965, -0.779088),
--    vpoint(0.388830, 0.381708, -0.651643),
--    vpoint(0.666557, 0.363844, -0.373362),
--    vpoint(0.777647, 0.311073, 0.010863),
--    vpoint(0.668031, 0.311074, 0.398231),
--    vpoint(0.376527, 0.374131, 0.661485),
--    vpoint(-0.010862, 0.311073, 0.777647),
--    vpoint(-0.390574, 0.363844, 0.656619),
--    vpoint(-0.678894, 0.311074, 0.379417),
--    vpoint(-0.777647, 0.311074, -0.010862),
--    vpoint(-0.656824, 0.374131, -0.384599),
--    vpoint(-0.379417, 0.311074, -0.678893),
--    vpoint(0.112758, 1.095153, 0.473684),
--    vpoint(0.189389, 1.097969, 0.375961),
--    vpoint(0.231078, 1.039332, 0.486629),
--    vpoint(0.212252, 1.143772, 0.593166),
--    vpoint(0.321897, 1.145405, 0.438386),
--    vpoint(0.253915, 1.236605, 0.519741),
--    vpoint(0.171172, 1.241927, 0.423361),
--    vpoint(-0.212785, 1.132481, 0.594174),
--    vpoint(-0.311047, 1.144340, 0.448341),
--    vpoint(-0.202055, 1.165019, 0.358402),
--    vpoint(-0.218891, 1.251477, 0.486628),
--    vpoint(-0.201075, 1.033039, 0.473684),
--    vpoint(-0.085836, 1.169609, 0.451970)}

-- rabbitNormals = {
--     vpoint(-0.2645, 0.8208, 0.5063),
--     vpoint(0.3750, 0.3463, -0.8599),
--     vpoint(0.1924, -0.6113, -0.7677),
--     vpoint(0.1488, 0.1210, -0.9814),
--     vpoint(0.3887, -0.8392, -0.3803),
--     vpoint(0.8856, -0.1074, -0.4519),
--     vpoint(0.5229, -0.8301, 0.1935),
--     vpoint(0.0398, -0.8273, 0.5604),
--     vpoint(0.1560, -0.9330, 0.3243),
--     vpoint(0.2225, 0.2616, 0.9392),
--     vpoint(0.0889, 0.9245, 0.3706),
--     vpoint(-0.1047, -0.9943, 0.0178),
--     vpoint(-0.3349, 0.9005, 0.2774),
--     vpoint(-0.7416, -0.2517, 0.6219),
--     vpoint(-0.5168, -0.8217, 0.2405),
--     vpoint(-0.5378, 0.8233, -0.1814),
--     vpoint(-0.8065, 0.5596, -0.1907),
--     vpoint(-0.7362, 0.2755, -0.6182),
--     vpoint(-0.0856, -0.9894, -0.1177),
--     vpoint(-0.7051, -0.5969, -0.3828),
--     vpoint(-0.3608, -0.1719, -0.9166),
--     vpoint(-0.2588, 0.0117, -0.9659),
--     vpoint(-0.3445, 0.4289, -0.8351),
--     vpoint(0.4067, 0.4508, -0.7946),
--     vpoint(0.6958, 0.2251, -0.6820),
--     vpoint(0.6952, -0.1323, 0.7065),
--     vpoint(-0.6958, 0.2251, 0.6820),
--     vpoint(-0.6960, -0.1323, -0.7057),
--     vpoint(0.8551, 0.2255, 0.4669),
--     vpoint(-0.1647, 0.3035, 0.9385),
--     vpoint(-0.3443, 0.0464, -0.9377),
--     vpoint(-0.6867, -0.1843, -0.7032),
--     vpoint(-0.6645, 0.7419, -0.0897),
--     vpoint(-0.1282, 0.6966, 0.7060),
--     vpoint(0.5306, 0.7419, 0.4099),
--     vpoint(0.1004, 0.6835, -0.7230),
--     vpoint(0.7492, -0.4020, 0.5265),
--     vpoint(-0.3755, -0.8700, -0.3195),
--     vpoint(0.5212, 0.7507, -0.4060),
--     vpoint(-0.5552, 0.6548, 0.5128),
--     vpoint(-0.7656, 0.4592, 0.4506),
--     vpoint(0.8121, 0.0913, 0.5763),
--     vpoint(0.3019, 0.3324, -0.8935),
--     vpoint(-0.7842, 0.1023, -0.6120),
--     vpoint(0.4982, -0.6744, -0.5450),
--     vpoint(-0.3923, -0.6906, 0.6077),
--     vpoint(0.7021, -0.5142, 0.4926),
--     vpoint(0.5221, 0.5864, 0.6192),
--     vpoint(0.3616, 0.7943, -0.4881),
--     vpoint(-0.5500, 0.6572, -0.5154),
--     vpoint(-0.7284, 0.4384, 0.5265),
--     vpoint(0.5355, -0.5533, -0.6380),
--     vpoint(-0.4435, -0.5915, -0.6734),
--     vpoint(-0.6919, -0.5884, 0.4184),
--     vpoint(0.7070, 0.0201, -0.7070),
--     vpoint(0.9658, 0.0131, -0.2588),
--     vpoint(0.9621, 0.0130, 0.2723),
--     vpoint(0.6726, 0.0208, 0.7397),
--     vpoint(0.2588, 0.0131, 0.9658),
--     vpoint(-0.3011, 0.0200, 0.9534),
--     vpoint(-0.7070, 0.0131, 0.7070),
--     vpoint(-0.9694, 0.0130, 0.2453),
--     vpoint(-0.9523, 0.0208, -0.3044),
--     vpoint(0.0091, 0.9999, 0.0056),
--     vpoint(-0.7070, 0.0131, -0.7070),
--     vpoint(0.3124, 0.0228, -0.9497),
--     vpoint(0.0000, -1.0000, -0.0000),
--     vpoint(0.2961, 0.8914, -0.3433),
--     vpoint(0.6368, 0.5171, 0.5719),
--     vpoint(-0.4453, 0.8709, -0.2078),
--     vpoint(0.3880, 0.8258, 0.4093),
--     vpoint(0.3368, 0.6405, -0.6901),
--     vpoint(0.3559, -0.3704, -0.8580),
--     vpoint(0.6291, 0.7365, -0.2485),
--     vpoint(0.8243, -0.1001, -0.5573),
--     vpoint(0.6411, -0.7537, -0.1443),
--     vpoint(0.8183, 0.5659, -0.1011),
--     vpoint(0.9190, -0.3348, -0.2082),
--     vpoint(0.2299, 0.9732, -0.0083),
--     vpoint(0.9361, -0.1308, 0.3265),
--     vpoint(0.8361, 0.4778, 0.2697),
--     vpoint(0.6393, -0.3524, 0.6835),
--     vpoint(0.2761, -0.3897, 0.8786),
--     vpoint(0.7527, -0.2187, 0.6210),
--     vpoint(0.1211, 0.8950, 0.4294),
--     vpoint(-0.2037, 0.3969, 0.8950),
--     vpoint(-0.2835, -0.2843, 0.9159),
--     vpoint(-0.6184, 0.4608, 0.6366),
--     vpoint(-0.6209, -0.3817, 0.6847),
--     vpoint(-0.1839, 0.8629, -0.4707),
--     vpoint(-0.4969, -0.7343, 0.4625),
--     vpoint(-0.9417, -0.1518, 0.3003),
--     vpoint(-0.8644, 0.4002, 0.3045),
--     vpoint(-0.8483, 0.3440, -0.4025),
--     vpoint(-0.7808, -0.5733, -0.2484),
--     vpoint(-0.9552, -0.2859, -0.0760),
--     vpoint(-0.3915, -0.7708, -0.5026),
--     vpoint(-0.5885, -0.1674, -0.7910),
--     vpoint(-0.1378, 0.1267, -0.9823),
--     vpoint(0.0867, -0.9954, -0.0416),
--     vpoint(-0.0693, -0.6202, -0.7814),
--     vpoint(0.6424, 0.5752, 0.5064),
--     vpoint(-0.6424, 0.5752, -0.5064),
--     vpoint(0.0359, -0.9951, -0.0918),
--     vpoint(0.0628, 0.0346, -0.9974),
--     vpoint(0.3086, -0.1565, -0.9382),
--     vpoint(0.9430, -0.0711, -0.3251),
--     vpoint(0.2231, 0.0769, 0.9718),
--     vpoint(-0.0634, 0.0347, 0.9974),
--     vpoint(-0.3076, -0.1565, 0.9386),
--     vpoint(-0.9431, -0.0712, 0.3249),
--     vpoint(-0.2226, 0.0764, -0.9719),
--     vpoint(-0.3802, 0.5428, 0.7489),
--     vpoint(0.5994, -0.1807, 0.7798),
--     vpoint(0.0537, 0.0591, -0.9968),
--     vpoint(0.3151, -0.0885, -0.9449),
--     vpoint(0.9130, -0.0857, -0.3988),
--     vpoint(0.8470, 0.2927, 0.4437),
--     vpoint(-0.2642, -0.0895, 0.9603),
--     vpoint(-0.8988, -0.0799, 0.4309),
--     vpoint(-0.8862, 0.2771, -0.3713),
--     vpoint(-0.1892, 0.5986, -0.7784),
--     vpoint(-0.4477, 0.6892, -0.5697),
--     vpoint(-0.7568, 0.6250, 0.1915),
--     vpoint(-0.5792, 0.6614, 0.4765),
--     vpoint(0.1115, 0.6688, 0.7350),
--     vpoint(0.7514, 0.6246, 0.2126),
--     vpoint(0.7401, 0.6611, -0.1232),
--     vpoint(-0.2507, 0.0097, -0.9680),
--     vpoint(0.7083, 0.0196, -0.7056),
--     vpoint(0.9611, 0.0173, -0.2755),
--     vpoint(0.9659, 0.0098, 0.2588),
--     vpoint(0.7071, 0.0098, 0.7071),
--     vpoint(0.2842, 0.0196, 0.9586),
--     vpoint(-0.2588, 0.0098, 0.9659),
--     vpoint(-0.6946, 0.0173, 0.7192),
--     vpoint(-0.9659, 0.0098, 0.2588),
--     vpoint(-0.9659, 0.0098, -0.2588),
--     vpoint(-0.0048, 1.0000, 0.0083),
--     vpoint(-0.7254, 0.0196, -0.6880),
--     vpoint(0.2588, 0.0096, -0.9659),
--     vpoint(0.4034, 0.8462, -0.3482),
--     vpoint(0.6008, 0.3956, 0.6946),
--     vpoint(0.1921, -0.4374, 0.8785),
--     vpoint(-0.1909, 0.3829, 0.9039),
--     vpoint(-0.6007, 0.4396, 0.6678),
--     vpoint(-0.1190, 0.7859, -0.6068),
--     vpoint(-0.0830, -0.9945, 0.0632),
--     vpoint(-0.0481, -0.9980, -0.0404),
--     vpoint(0.4603, -0.0894, 0.8832),
--     vpoint(-0.1402, -0.0278, 0.9897)}

-- rabbitTriangles = {
--     {{4,1,1}, {5,2,1}, {75,3,1}, color_white},
--     {{17,4,2}, {22,5,2}, {19,6,2}, color_white},
--     {{18,7,3}, {23,8,3}, {20,9,3}, color_white},
--     {{17,4,4}, {19,6,4}, {18,7,4}, color_white},
--     {{20,9,5}, {23,8,5}, {31,10,5}, color_white},
--     {{22,5,6}, {25,11,6}, {23,8,6}, color_white},
--     {{26,12,7}, {30,13,7}, {31,10,7}, color_white},
--     {{30,13,8}, {34,14,8}, {37,15,8}, color_white},
--     {{30,13,9}, {37,15,9}, {31,16,9}, color_white},
--     {{29,17,10}, {32,18,10}, {33,19,10}, color_white},
--     {{21,20,11}, {35,21,11}, {28,22,11}, color_white},
--     {{31,23,12}, {37,15,12}, {44,24,12}, color_white},
--     {{35,21,13}, {21,20,13}, {39,25,13}, color_white},
--     {{36,26,14}, {38,27,14}, {41,28,14}, color_white},
--     {{41,28,15}, {44,24,15}, {37,15,15}, color_white},
--     {{39,29,16}, {21,30,16}, {42,31,16}, color_white},
--     {{39,29,17}, {42,31,17}, {40,32,17}, color_white},
--     {{42,31,18}, {45,33,18}, {43,34,18}, color_white},
--     {{31,35,19}, {44,36,19}, {20,9,19}, color_white},
--     {{43,34,20}, {46,37,20}, {44,36,20}, color_white},
--     {{45,33,21}, {18,7,21}, {46,37,21}, color_white},
--     {{16,38,22}, {68,39,22}, {1,40,22}, color_white},
--     {{42,31,23}, {17,4,23}, {45,33,23}, color_white},
--     {{49,41,24}, {47,42,24}, {60,43,24}, color_white},
--     {{54,44,25}, {53,45,25}, {48,46,25}, color_white},
--     {{55,47,26}, {52,48,26}, {54,44,26}, color_white},
--     {{58,49,27}, {56,50,27}, {50,51,27}, color_white},
--     {{51,52,28}, {57,53,28}, {58,49,28}, color_white},
--     {{47,42,29}, {64,54,29}, {60,43,29}, color_white},
--     {{64,54,30}, {47,42,30}, {65,55,30}, color_white},
--     {{49,56,31}, {59,57,31}, {67,58,31}, color_white},
--     {{59,57,32}, {66,59,32}, {67,58,32}, color_white},
--     {{78,60,33}, {77,61,33}, {5,2,33}, color_white},
--     {{75,3,34}, {74,62,34}, {4,1,34}, color_white},
--     {{73,63,35}, {72,64,35}, {3,65,35}, color_white},
--     {{2,66,36}, {69,67,36}, {68,68,36}, color_white},
--     {{83,69,37}, {82,70,37}, {84,71,37}, color_white},
--     {{82,72,38}, {80,73,38}, {81,74,38}, color_white},
--     {{85,75,39}, {84,71,39}, {86,76,39}, color_white},
--     {{86,77,40}, {83,69,40}, {85,78,40}, color_white},
--     {{80,79,41}, {83,80,41}, {86,76,41}, color_white},
--     {{83,69,42}, {84,71,42}, {85,78,42}, color_white},
--     {{84,71,43}, {81,81,43}, {86,76,43}, color_white},
--     {{81,81,44}, {80,79,44}, {86,76,44}, color_white},
--     {{82,82,45}, {81,81,45}, {84,71,45}, color_white},
--     {{82,83,46}, {83,80,46}, {80,79,46}, color_white},
--     {{87,84,47}, {91,85,47}, {92,86,47}, color_white},
--     {{87,84,48}, {92,86,48}, {90,87,48}, color_white},
--     {{92,88,49}, {89,89,49}, {90,90,49}, color_white},
--     {{89,89,50}, {88,91,50}, {90,92,50}, color_white},
--     {{88,91,51}, {87,93,51}, {90,94,51}, color_white},
--     {{91,95,52}, {89,89,52}, {92,86,52}, color_white},
--     {{91,96,53}, {88,91,53}, {89,89,53}, color_white},
--     {{91,97,54}, {87,93,54}, {88,91,54}, color_white},
--     {{6,98,55}, {70,99,55}, {7,100,55}, color_white},
--     {{7,100,56}, {71,101,56}, {8,102,56}, color_white},
--     {{9,103,57}, {71,101,57}, {72,64,57}, color_white},
--     {{10,104,58}, {72,64,58}, {73,63,58}, color_white},
--     {{10,104,59}, {74,62,59}, {11,105,59}, color_white},
--     {{12,106,60}, {74,62,60}, {75,3,60}, color_white},
--     {{12,106,61}, {76,107,61}, {13,108,61}, color_white},
--     {{14,109,62}, {76,107,62}, {77,61,62}, color_white},
--     {{15,110,63}, {77,61,63}, {78,60,63}, color_white},
--     {{4,111,64}, {2,112,64}, {5,113,64}, color_white},
--     {{15,110,65}, {79,114,65}, {16,38,65}, color_white},
--     {{6,98,66}, {68,68,66}, {69,67,66}, color_white},
--     {{11,115,67}, {15,116,67}, {7,117,67}, color_white},
--     {{3,65,68}, {69,67,68}, {2,66,68}, color_white},
--     {{32,18,69}, {27,118,69}, {28,22,69}, color_white},
--     {{5,2,70}, {2,119,70}, {78,60,70}, color_white},
--     {{3,65,71}, {4,1,71}, {73,63,71}, color_white},
--     {{17,4,72}, {24,120,72}, {22,5,72}, color_white},
--     {{18,7,73}, {19,6,73}, {23,8,73}, color_white},
--     {{24,120,74}, {28,22,74}, {22,5,74}, color_white},
--     {{19,6,75}, {22,5,75}, {23,8,75}, color_white},
--     {{31,10,76}, {23,8,76}, {26,12,76}, color_white},
--     {{22,5,77}, {28,22,77}, {25,11,77}, color_white},
--     {{23,8,78}, {25,11,78}, {26,12,78}, color_white},
--     {{24,120,79}, {21,20,79}, {28,22,79}, color_white},
--     {{25,11,80}, {27,118,80}, {26,12,80}, color_white},
--     {{25,11,81}, {28,22,81}, {27,118,81}, color_white},
--     {{26,12,82}, {29,17,82}, {30,13,82}, color_white},
--     {{33,19,83}, {30,13,83}, {29,17,83}, color_white},
--     {{26,12,84}, {27,118,84}, {29,17,84}, color_white},
--     {{28,22,85}, {35,21,85}, {32,18,85}, color_white},
--     {{36,26,86}, {32,18,86}, {35,21,86}, color_white},
--     {{33,19,87}, {36,26,87}, {34,14,87}, color_white},
--     {{38,27,88}, {35,21,88}, {39,25,88}, color_white},
--     {{34,14,89}, {36,26,89}, {41,28,89}, color_white},
--     {{24,120,90}, {42,31,90}, {21,30,90}, color_white},
--     {{34,14,91}, {41,28,91}, {37,15,91}, color_white},
--     {{41,28,92}, {38,27,92}, {40,121,92}, color_white},
--     {{38,27,93}, {39,25,93}, {40,121,93}, color_white},
--     {{40,32,94}, {42,31,94}, {43,34,94}, color_white},
--     {{41,122,95}, {43,34,95}, {44,36,95}, color_white},
--     {{41,122,96}, {40,32,96}, {43,34,96}, color_white},
--     {{44,36,97}, {46,37,97}, {20,9,97}, color_white},
--     {{43,34,98}, {45,33,98}, {46,37,98}, color_white},
--     {{45,33,99}, {17,4,99}, {18,7,99}, color_white},
--     {{51,52,100}, {55,47,100}, {57,53,100}, color_white},
--     {{46,37,101}, {18,7,101}, {20,9,101}, color_white},
--     {{48,46,102}, {50,51,102}, {54,44,102}, color_white},
--     {{50,51,103}, {48,123,103}, {58,49,103}, color_white},
--     {{59,57,104}, {63,124,104}, {66,59,104}, color_white},
--     {{48,46,105}, {53,45,105}, {51,125,105}, color_white},
--     {{51,125,106}, {53,45,106}, {52,48,106}, color_white},
--     {{53,45,107}, {54,44,107}, {52,48,107}, color_white},
--     {{54,44,108}, {50,51,108}, {55,47,108}, color_white},
--     {{50,51,109}, {56,50,109}, {55,47,109}, color_white},
--     {{55,47,110}, {56,50,110}, {57,53,110}, color_white},
--     {{56,50,111}, {58,49,111}, {57,53,111}, color_white},
--     {{58,49,112}, {48,123,112}, {51,52,112}, color_white},
--     {{65,55,113}, {47,42,113}, {49,56,113}, color_white},
--     {{62,126,114}, {63,124,114}, {61,127,114}, color_white},
--     {{49,41,115}, {60,43,115}, {59,128,115}, color_white},
--     {{59,128,116}, {60,43,116}, {61,127,116}, color_white},
--     {{60,43,117}, {62,126,117}, {61,127,117}, color_white},
--     {{60,43,118}, {64,54,118}, {62,126,118}, color_white},
--     {{63,124,119}, {65,55,119}, {66,59,119}, color_white},
--     {{65,55,120}, {67,58,120}, {66,59,120}, color_white},
--     {{65,55,121}, {49,56,121}, {67,58,121}, color_white},
--     {{2,119,122}, {68,39,122}, {79,114,122}, color_white},
--     {{2,119,123}, {79,114,123}, {78,60,123}, color_white},
--     {{5,2,124}, {77,61,124}, {76,107,124}, color_white},
--     {{5,2,125}, {76,107,125}, {75,3,125}, color_white},
--     {{4,1,126}, {74,62,126}, {73,63,126}, color_white},
--     {{3,65,127}, {72,64,127}, {71,101,127}, color_white},
--     {{3,65,128}, {71,101,128}, {70,99,128}, color_white},
--     {{16,38,129}, {79,114,129}, {68,39,129}, color_white},
--     {{6,98,130}, {69,67,130}, {70,99,130}, color_white},
--     {{7,100,131}, {70,99,131}, {71,101,131}, color_white},
--     {{9,103,132}, {8,102,132}, {71,101,132}, color_white},
--     {{10,104,133}, {9,103,133}, {72,64,133}, color_white},
--     {{10,104,134}, {73,63,134}, {74,62,134}, color_white},
--     {{12,106,135}, {11,105,135}, {74,62,135}, color_white},
--     {{12,106,136}, {75,3,136}, {76,107,136}, color_white},
--     {{14,109,137}, {13,108,137}, {76,107,137}, color_white},
--     {{15,110,138}, {14,109,138}, {77,61,138}, color_white},
--     {{4,111,139}, {3,129,139}, {2,112,139}, color_white},
--     {{15,110,140}, {78,60,140}, {79,114,140}, color_white},
--     {{6,98,141}, {1,130,141}, {68,68,141}, color_white},
--     {{7,117,67}, {8,131,67}, {9,132,67}, color_white},
--     {{9,132,67}, {10,133,67}, {7,117,67}, color_white},
--     {{10,133,67}, {11,115,67}, {7,117,67}, color_white},
--     {{11,115,67}, {12,134,67}, {13,135,67}, color_white},
--     {{13,135,67}, {14,136,67}, {15,116,67}, color_white},
--     {{15,116,67}, {16,137,67}, {1,138,67}, color_white},
--     {{1,138,67}, {6,139,67}, {15,116,67}, color_white},
--     {{6,139,67}, {7,117,67}, {15,116,67}, color_white},
--     {{11,115,67}, {13,135,67}, {15,116,67}, color_white},
--     {{3,65,142}, {70,99,142}, {69,67,142}, color_white},
--     {{32,18,143}, {29,17,143}, {27,118,143}, color_white},
--     {{33,19,144}, {34,14,144}, {30,13,144}, color_white},
--     {{36,26,145}, {33,19,145}, {32,18,145}, color_white},
--     {{38,27,146}, {36,26,146}, {35,21,146}, color_white},
--     {{24,120,147}, {17,4,147}, {42,31,147}, color_white},
--     {{51,52,148}, {52,48,148}, {55,47,148}, color_white},
--     {{59,57,149}, {61,127,149}, {63,124,149}, color_white},
--     {{62,126,150}, {64,54,150}, {63,124,150}, color_white},
--     {{63,124,151}, {64,54,151}, {65,55,151}, color_white}}

-- function rabbitMesh()
--    return meshdata(rabbitVertices, rabbitNormals, rabbitTriangles,
--       -- backface
--       true)
-- end
-------------------------------------------------------------------------------
-->8
-- borrowed

-- trifill by @p01
function p01_trapeze_h(l,r,lt,rt,y0,y1)
   lt,rt=(lt-l)/(y1-y0),(rt-r)/(y1-y0)
   if(y0<0)l,r,y0=l-y0*lt,r-y0*rt,0 
   for y0=y0,min(y1,128) do
      rectfill(l,y0,r,y0)
      l+=lt
      r+=rt
   end
   end
function p01_trapeze_w(t,b,tt,bt,x0,x1)
   tt,bt=(tt-t)/(x1-x0),(bt-b)/(x1-x0)
   if(x0<0)t,b,x0=t-x0*tt,b-x0*bt,0 
   for x0=x0,min(x1,128) do
      rectfill(x0,t,x0,b)
      t+=tt
      b+=bt
   end
   end

function trifill(x0,y0,x1,y1,x2,y2,c)
   color(c)
   if(y1<y0)x0,x1,y0,y1=x1,x0,y1,y0
   if(y2<y0)x0,x2,y0,y2=x2,x0,y2,y0
   if(y2<y1)x1,x2,y1,y2=x2,x1,y2,y1
   if max(x2,max(x1,x0))-min(x2,min(x1,x0)) > y2-y0 then
      local col=x0+(x2-x0)/(y2-y0)*(y1-y0)
      p01_trapeze_h(x0,x0,x1,col,y0,y1)
      p01_trapeze_h(x1,col,x2,x2,y1,y2)
   else
      if(x1<x0)x0,x1,y0,y1=x1,x0,y1,y0
      if(x2<x0)x0,x2,y0,y2=x2,x0,y2,y0
      if(x2<x1)x1,x2,y1,y2=x2,x1,y2,y1
      local col=y0+(y2-y0)/(x2-x0)*(x1-x0)
      p01_trapeze_w(y0,y0,y1,col,x0,x1)
      p01_trapeze_w(y1,col,y2,y2,x1,x2)
      end
   end

--quicksort from three-d by jimmi
function gt(_x,_y) return (_x > _y) end
function lt(_x,_y) return (_x < _y) end
function qsort(_data,_lo,_hi,_key,_cmp)
  _cmp = _cmp or gt
  if _lo < _hi then
    local pivot, p = _data[_hi][_key], _lo
    for j=_lo,_hi-1 do
      if _cmp(_data[j][_key], pivot) then
        _data[p], _data[j], p = _data[j], _data[p], p+1
      end
    end
    _data[p], _data[_hi] = _data[_hi], _data[p]
    qsort(_data,_lo,p-1,_key,_cmp)
    qsort(_data,p+1,_hi,_key,_cmp)
  end
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

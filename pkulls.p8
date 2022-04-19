pico-8 cartridge // http://www.pico-8.com
version 35
__lua__

-- main
dt=1/30
g=-9.8/30 -- gravity for 1 frame
p={}
boats={}
wind={}
waves={}
windtrails={}
cannonballs={}
next_windrail=0
was_colliding=false
col_grace_period=0
camshake=nil
metervalue=0
boarding_state=0
boarding_state_t=0

should_debug_col = false
debug_coll={}

speed_mult=0.1
max_speed=2
rot_speed=0.01
max_cannonball_dist=40
max_cannonball_halfdist=max_cannonball_dist * 0.5
max_wind=5
c2c_meter_inc=0.002
c2c_dmg_inc=0.001
wave_sprites={4,5,6,5}
wave_sprites_count=count(wave_sprites)
splash_duration=45
splash_frames={10,11,12,13,14,15}
expl_duration=20
expl_frames={16,17,18,19,20,21,22}
c2c_frames={23,23,24,24,25,25,26,26,27,27}
tethering_frames={64,68}
boarding_frames={72}
loot_frames={128,130,132,134}

function _init()
   p=make_boat(64, 64)
   wind={dir=rnd(1), str=2+rnd(max_wind-2)}
   for i=1,30 do
      add(waves, make_wave(64,64))
   end
   for i=1,20 do
      add(boats, make_boat(rnd(256)-128, rnd(256)-128))
   end

   pal(2, 132, 1) -- swap purple for dark brown
   pal(3, 140, 1) -- swap dark green for dark blue
   pal(14, 128, 1) -- swap pink for very dark brown
end

function _update()
   update_boarding()
   
   if boarding_state == 0 then
      update_wind()
      update_inputs()
      update_boat(p)
      update_cannonballs()
      update_waves()
      update_windtrails()
      update_camshake()
   end
end

function _draw()
   cls(3)

   camdx, camdy = getcamshake()
   camera(p.x-64+camdx, p.y-64+camdy)
   draw_waves()
   draw_cannonballs_splash()
   draw_ai_boats()
   draw_player_boat(p)
   draw_cannonballs()
   draw_windtrails()
   draw_debug_col()

   camera()
   draw_hud()
   draw_boarding()

   -- print("wind: "..wind.dir.." "..wind.str)
   -- print("p: "..p.speed.." "..p.state)
   -- print("bstate: "..boarding_state.." t: "..boarding_state_t)
   -- if count(cannonballs) > 0 then
   --    local ball = cannonballs[1]
   --    print("[ball] x: "..ball.x.." y: "..ball.y.." d: "..ball.dist.." state: "..ball.state)
   -- end
end

-- inputs and controls

-- state 0: ok 1: sinked
function make_boat(x,y)
   local scale=1.5
   return {x=x, y=y, speed=0, v=0, dir=0, rx=0, ry=0, w=8*scale, h=8*scale, scale=scale, state=0, hp=1}
end
function update_inputs()
   if btn(0) then -- left
      p.dir += rot_speed
   elseif btn(1) then -- right
      p.dir -= rot_speed
   end
   if btnp(2) then -- up
      p.speed += 1
   elseif btnp(3) then -- down
      p.speed -= 1
   end
   if (btnp(4) and btn(5)) or (btn(4) and btnp(5)) then
      try_boarding()
   elseif not was_colliding then -- no cannons while colliding, too close!
      if btnp(4) then -- square
         make_cannonball(p, 1)
      end
      if btnp(5) then -- X
         make_cannonball(p, -1)
      end
   end
   p.dir = normalize_angle(p.dir)
   p.speed = min(max_speed, max(0, p.speed))
end

function update_boat(boat)
   local delta_angle = shortest_angle_between_normalized(boat.dir, wind.dir)
   -- boat movement vs wind inspired from:
   -- https://sailing-blog.nauticed.org/sailboat-speed-versus-sailing-angle/
   local boat_min_speed = 0.5 * wind.str * boat.speed * speed_mult
   local boat_max_speed = 1.5 * wind.str * boat.speed * speed_mult
   boat.v = boat_min_speed + (boat_max_speed - boat_min_speed) * -1 * sin(delta_angle)
   -- apply only integer delta values and re-apply the rest the next frame
   local dx, dy = p.rx + cos(boat.dir) * boat.v, p.ry + sin(boat.dir) * boat.v
   local fdx, fdy = flr(dx), flr(dy)
   local col = check_collisions(boat.x+fdx, boat.y+fdy, 0, 0, boat.w, boat.h, 0.75)
   if col then
      if not was_colliding and col_grace_period == 0 then
         add_camshake(8,0,1,0,0.3)
      end
      if boat == p then
         metervalue = min(1, metervalue + c2c_meter_inc)
         apply_damage(p, c2c_dmg_inc, col)
         apply_damage(col, c2c_dmg_inc, p)
      end
      was_colliding = col
   else
      p.rx, p.ry = dx - fdx, dy - fdy
      boat.x += fdx
      boat.y += fdy
      if was_colliding then
         was_colliding = false
         col_grace_period += 15
      end
   end
   col_grace_period = max(0, col_grace_period-1)
end

boat_sprite_map = {x=0, y=0}
sail_sprites_map = {{x=2,y=0},{x=4,y=0},{x=6,y=0}}
function draw_boat(boat)
   if boat.state == 0 then -- state ok
      draw_rotated_tile(boat.x, boat.y, -boat.dir, boat_sprite_map.x, boat_sprite_map.y, 1, false, boat.scale)
      local sail_sprite_map = sail_sprites_map[boat.speed+1]
      draw_rotated_tile(boat.x, boat.y, -wind.dir, sail_sprite_map.x, sail_sprite_map.y, 1, false, boat.scale)
   elseif boat.state == 1 then -- state sunk
      scaled_spr(2, boat.x, boat.y, boat.scale)
   end
end

c2c_frame=1
function draw_player_boat()
   draw_boat(p)

   if was_colliding and boarding_state == 0 then
      scaled_spr(c2c_frames[c2c_frame], p.x, p.y, p.scale)
      scaled_spr(c2c_frames[c2c_frame], was_colliding.x, was_colliding.y, p.scale)
      c2c_frame = (c2c_frame+1) % count(c2c_frames) + 1
   end
end

function draw_hud()
   hpx,hpy = 2,2
   w,h = 20,2
   local hpcol = p.hp > 0.5 and 11 or 8
   rectfill(hpx, hpy, hpx+w, hpy+h, 1)
   rectfill(hpx, hpy, hpx+w*p.hp,hpy+h, hpcol)

   meterx,metery = 2,4
   local ismeter_full = metervalue >= 1
   metercol = ismeter_full and 9 or 10
   rectfill(meterx, metery, meterx+w, metery+h, 1)
   if metervalue > 0 then
      rectfill(meterx, metery, meterx+w*metervalue,metery+h, metercol)
   end
   if ismeter_full and was_colliding then
      print("❎🅾️ to board", 40, 122)
   end
end

function draw_ai_boats()
   for boat in all(boats) do
      pal(4,2,0)
      draw_boat(boat)
      
      if boat.state == 0 and boat.hp < 1 then
         local offsetx, offsety = 4*boat.scale, 3*boat.scale
         local barx, bary = boat.x - offsetx, boat.y + offsety
         local w,h = 8*boat.scale, 1*boat.scale
         rectfill(barx, bary, barx + w, bary + h, 1)
         if boat.hp > 0 then
            local hpcol = boat.hp > 0.5 and 11 or 8
            rectfill(barx, bary, barx + w*boat.hp, bary + h, hpcol)
         end
      end
      pal(0)
   end
end

function check_collisions(x, y, dx, dy, w, h, s)
   s = s or 1
   for b in all(boats) do
      hw,hh, bhw, bhh = s*w*0.5, s*h*0.5, s*b.w*0.5, s*b.h*0.5
      local res = intersect(x-hw,y-hh,x+dx+hw,y+dy+hh, b.x-bhw,b.y-bhh,b.x+bhw,b.y+bhh)
      if should_debug_col then
         add(debug_coll, {x-hw,y-hh,x+dx+hw,y+dy+hh})
         add(debug_coll, {b.x-bhw,b.y-bhh,b.x+bhw,b.y+bhh})
      end
      if res and b.state != 1 then -- boat not sunk
         return b
      end
   end
   return false
end

function add_camshake(freqx, freqy, ampx, ampy, duration)
   camshake={dx=0, dy=0, t=0, freqx=freqx/30, freqy=freqy/30, ampx=ampx, ampy=ampy, lifetime=flr(duration*30)}
end

function getcamshake()
   if (camshake == nil) return 0,0
   return camshake.dx, camshake.dy
end

function update_camshake()
   if (camshake == nil) return

   if camshake.t >= camshake.lifetime then
      camshake = nil
      return
   end
   camshake.dx = -camshake.ampx * sin(camshake.t * camshake.freqx)
   camshake.dy = -camshake.ampy * sin(camshake.t * camshake.freqy)
   camshake.t += 1
end

function apply_damage(boat, dmg, instigator)
   boat.hp = max(0, boat.hp-dmg)
   if boat.hp == 0 then
      boat.state = 1 -- sink boat
   end
end

function try_boarding()
   if (metervalue < 1) return
   if (not was_colliding) return
   local target = was_colliding -- last collision is the target
   set_boarding_state(1) -- tethering
   metervalue = 0
end
function set_boarding_state(state)
   if boarding_state == 3 and state == 0 then
      apply_damage(was_colliding, 1.0) -- on exit loot, kill target
   end

   boarding_state=state
   boarding_state_t = 0
end

function update_boarding()
   if (boarding_state == 0) return -- no boarding
   if boarding_state == 1 then
      if (boarding_state_t == 120) set_boarding_state(2) -- to boarding outcome
   elseif boarding_state == 2 then
      if (boarding_state_t == 120) set_boarding_state(3) -- to loot
   elseif boarding_state == 3 then
      if (boarding_state_t == 120) set_boarding_state(0) -- to no boarding
   end
   boarding_state_t += 1
end

function draw_boarding()
   if boarding_state == 1 then
      local spriteidx = (boarding_state_t \ 15) % 2 + 1
      local spritenum = tethering_frames[spriteidx]
      scaled_spr(spritenum, 64, 64, 4, 4, 3)
      rect(16, 16, 112, 112, 10)
      rect(15, 15, 113, 113, 10)
      if boarding_state_t > 30 then
         rect(14, 14, 114, 114, 9)
         rect(13, 13, 115, 115, 9)
         if boarding_state_t > 60 then
            rect(12, 12, 116, 116, 8)
            rect(11, 11, 117, 117, 8)
            if boarding_state_t > 90 then
               rect(10, 10, 118, 118, 4)
               rect(9, 9, 119, 119, 4)
            end
         end
      end
   elseif boarding_state == 2 then -- outcome
      rectfill(9, 9, 119, 119, 8)
      scaled_spr(boarding_frames[1], 64, 64, 4, 4, 3)
   elseif boarding_state == 3 then -- loot
      rectfill(9, 9, 119, 119, 8)
      local spriteidx = flr(boarding_state_t / 121 * count(loot_frames)) + 1
      local spritenum = loot_frames[spriteidx]
      scaled_spr(spritenum, 64, 64, 2, 2, 6)
   end
end

-- states: [0: inair] [1: splash] [2: hit]
function make_cannonball(boat, side)
   local dir = normalize_angle(boat.dir+side*0.25)
   add(cannonballs, {x=boat.x, y=boat.y, dir=dir, v=1.0, dist=0, state=0, state_t=0, r=2, instigator=boat})
end

function update_cannonballs()
   local toremove={}
   for i=1,count(cannonballs) do
      local ball = cannonballs[i]
      if ball.state == 0 then -- inair
         local dx, dy = ball.v*cos(ball.dir), ball.v*sin(ball.dir)
         local collision_res = check_collisions(ball.x, ball.y, dx, dy, ball.r, ball.r)
         if collision_res and invlerp(ball.dist, 0, max_cannonball_dist) > 0.8 then
            ball.x, ball.y = collision_res.x, collision_res.y
            ball.state = 2 -- -> hit
            apply_damage(collision_res, 0.2, ball.instigator)
            if ball.instigator == p then
               metervalue = min(1, metervalue+0.1)
            end
         else
            ball.x += dx
            ball.y += dy
            ball.dist += sqrt(dx*dx + dy*dy)
            if ball.dist > max_cannonball_dist then
               ball.state = 1 -- -> splash
            end
         end
      elseif ball.state == 1 or ball.state == 2 then -- splash or hit
         local durations = {splash_duration, expl_duration}
         local duration = durations[ball.state]
         if ball.state_t < duration then
            ball.state_t +=1
         else
            add(toremove, i)
         end
      end
   end
   for i in all(toremove) do
      deli(cannonballs, i)
   end
end

function draw_cannonballs()
   for ball in all(cannonballs) do
      if ball.state == 0 then -- inair
         local middist = abs(ball.dist - max_cannonball_halfdist)
         local r = lerp(ball.r,ball.r*2.5,invlerp(middist,max_cannonball_halfdist,0))
         circfill(ball.x, ball.y, r, 5)
      elseif ball.state == 2 then -- hit
         local expl_sprite_idx = flr(lerp(1, count(expl_frames)+0.5, invlerp(ball.state_t, 0, expl_duration)))
         local expl_sprite = expl_frames[expl_sprite_idx]
         scaled_spr(expl_sprite, ball.x, ball.y, 1, 1, 3)
      end
   end
end

function draw_cannonballs_splash()
   for ball in all(cannonballs) do
      if ball.state == 1 then -- splash
         local splash_sprite_idx = flr(lerp(1, count(splash_frames)+0.5, invlerp(ball.state_t, 0, splash_duration)))
         local splash_sprite = splash_frames[splash_sprite_idx]
         scaled_spr(splash_sprite, ball.x, ball.y, 1, 1, 2)
      end
   end
end


-- water and wind effects

function make_wave(cx, cy)
   local speed = rnd(90)+60
   return {x=cx+rnd(128)-64, y=cy+rnd(128)-64, state=0, substate=rnd(speed), speed=speed}
end


function update_waves()
   for i=1,count(waves) do
      w = waves[i]
      w.substate += 1
      if w.substate >= w.speed then
         w.substate = 0
         w.state += 1
      end
      -- loop the waves around the player
      local dx, dy = w.x - p.x, w.y - p.y
      if (dx > 64) w.x -= 128
      if (dx < -64) w.x += 128
      if (dy > 64) w.y -= 128
      if (dy < -64) w.y += 128
   end
end

function draw_waves()
   for w in all(waves) do
      spr(wave_sprites[w.state % wave_sprites_count +1], w.x, w.y)
   end
end

function update_wind()
   -- wind.dir = wind.dir + 0.001
   -- if wind.dir > 1 then
   --    wind.dir = 0
   -- end
end

function make_windtrail()
   add(windtrails, {x0=p.x+rnd(48)-24, y0=p.y+rnd(48)-24, life=0, length=0, maxlen=rnd(10)+3, maxdist=rnd(30)+15, dir=wind.dir})
end

function update_windtrails()
   next_windrail -= 1
   if next_windrail <=0 then
      make_windtrail()
      next_windrail = rnd(30)+30
   end

   local toremove={}
   local trailcount = count(windtrails)
   for i=1,trailcount do
      local w = windtrails[i]
      if w.length >= w.maxdist then
         add(toremove, i)
      else
         w.life +=1
         local mod = (max_wind - ceil(wind.str)) + 1
         if w.life % mod == 0 then
            w.length += 1
         end
      end
   end
   for i in all(toremove) do
      deli(windtrails, i)
   end
end

function draw_windtrails()
   for w in all(windtrails) do
      local l0 = max(0, w.length - w.maxlen)
      local x0,y0 = w.x0 + cos(w.dir) * l0, w.y0 + sin(w.dir) * l0
      local x1,y1 = w.x0 + cos(w.dir) * w.length, w.y0 + sin(w.dir) * w.length
      line(x0, y0, x1, y1, 6)
   end
end

function draw_debug_col()
   if should_debug_col then
      for c in all(debug_coll) do
         rect(c[1], c[2], c[3], c[4], 8)
      end
      debug_coll={}
   end
end


-- utils
function normalize_angle(a)
   local extra = flr(abs(a))
   if a > 1 then
      a -= extra
   elseif a < 0 then
      a += extra+1
   end
   return a
end

function shortest_angle_between_normalized(a1, a2)
   if a1 > a2 then
      a1,a2 = a2, a1
   end
   local delta = a2 - a1
   if delta < 0.5 then
      return delta
   else
      return 1.0 - delta
   end
end

function lerp(a, b, ratio)
   local delta = b-a
   return a + delta * ratio
end

function invlerp(val, a, b)
   return (val - a) / (b - a)
end

-- intersect rect a with rect b
function intersect(a_minx,a_miny,a_maxx,a_maxy, b_minx,b_miny,b_maxx,b_maxy)
   local xintersect = a_maxx > b_minx and a_minx < b_maxx
   if (not xintersect) return false
   local yintersect = a_maxy > b_miny and a_miny < b_maxy
   return yintersect

   -- if (not yintersect) return false
   -- local col_x, col_y = 0,0
   -- a_cx, a_cy = lerp(a_minx, a_maxx, 0.5), lerp(a_miny, a_maxy, 0.5)
   -- b_cx, b_cy = lerp(b_minx, b_maxx, 0.5), lerp(b_miny, b_maxy, 0.5)
   -- if a_cx < b_cx then
   --    col_x = b_minx - (a_cx - a_minx)
   -- else
   --    col_x = b_maxx + (a_cx - a_minx)
   -- end
   -- if a_cy < b_cy then
   --    col_y = b_miny - (a_cy - a_miny)
   -- else
   --    col_y = b_maxy + (a_cy - a_miny)
   -- end
   -- return {x= col_x, y= col_y}
end

-- @TheRoboZ https://www.lexaloffle.com/bbs/?pid=78451
function draw_rotated_tile(x,y,rot,mx,my,w,flip,scale)
  --print("draw_rotated_tile("..x..","..y..","..rot..","..mx..","..my..","..w..","..(flip and "true" or "false")..","..scale..")")
  scale = scale or 1
  w+=.8
  local halfw, cx  = scale*-w/2, mx + w/2 -.4
  local cs, ss, cy = cos(rot)/scale, -sin(rot)/scale, my-halfw/scale-.4
  local sx, sy, hx, hy = cx + cs*halfw, cy - ss*halfw, w*(flip and -4 or 4)*scale, w*4*scale

  --this just draw a bounding box to show the exact draw area
  --rect(x-hx,y-hy,x+hx,y+hy,5)

  for py = y-hy, y+hy do
    tline(x-hx, py, x+hx, py, sx + ss*halfw, sy + cs*halfw, cs/8, -ss/8)
    halfw+=1/8
  end
end

function scaled_spr(n, x, y, w, h, s)
   w,h,s = w or 1, h or 1, s or 1
   local x0,y0 = n*8, n*8
   local sx,sy = x0 % 128, y0 \ 128 * 8
   local sw,sh = 8*s*w, 8*s*h
   --printh("x0: "..x0.." y0: "..y0.." sx: "..sx.." sy: "..sy)
   sspr(sx, sy, w*8, h*8, x-sw*0.5, y-sh*0.5, sw, sh)
end


__gfx__
000000000000000000000470000000000000000000000000000000000000000000000000000f0000000000000000000000000000000cc00000077000000cc000
0000000000000000000047700000000000cc0000000000000000000000000000000000000007f0000000000000000000000cc00000c77c00007cc70000c33c00
007007004454544004004000000000000c00c00000cc000000cc000000000000000700000007700000000000000cc00000c77c000c7cc7c007c33c700c3333c0
00077000444444444444000000000000000000000c00c0000c00c000000700000004700000047000000cc00000c77c000c7cc7c0c7c33c7c7c3333c7c333333c
0007700044444444444004000000000000000cc000000cc000000000000700000004700000047000000cc00000c77c000c7cc7c0c7c33c7c7c3333c7c333333c
007007004454544004404404000000000000c00c0000c00c00000cc000000000000700000007700000000000000cc00000c77c000c7cc7c007c33c700c3333c0
0000000000000000cc004ccc0000000000000000000000000000c00c00000000000000000007f0000000000000000000000cc00000c77c00007cc70000c33c00
00000000000000000ccccc00000000000000000000000000000000000000000000000000000f0000000000000000000000000000000cc00000077000000cc000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000500000000000000050000000000000000000000000000000000
000000000077700000aaa000000a0000000900000004000000000000000000000000500000000000005050000050000000000000000000000000000000000000
00070000007a700000a9a00000a9a000009890000048400000040000000050000050050000005000000005000000000000000000000000000000000000000000
000000000077700000aaa000000a0000000900000004000000000000000000000000000000050000000000000000500000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
33333333355333333333333333333333333333333333355333333333333333334e44444e44444444444e44444444444e4e44444e444444444444444444444444
33333333335555333333333333333333333333333333335533333333333333334e44444e4444544444e4444444444ee44e44444e444454444444444444444444
3333333333333555333cc33333333333333333333333333553333333333333334e44444e44454444ee444444444ee4444e44444e444544444444444444444444
33333333333333355533ccccc3333333533333333333333335533333333333334e44444e4445444e4477777744e444444e44444e444544444477777744444444
53333333333cc333355533333333333355533333333333333c555333333333334e44444e444554e447777774ee4444444e44444e444554444777777444444444
355333cc3333cccc333553333333333333555333333333333ccc55c333333333ee44444eee445e444777777e44444444ee44444eee4454444777777444444444
3335553ccc33333333335555333333ee333355553cc3333333333553333333ee4e44444e4e4e54444fcfcfff4444444e4e44444e4e4454444fcfcfff44444444
33333355333333333333333553313e443333333555cccc333333335533333e44444444444ee455544ffffff444444ee4444444444e4455544ffffff444444444
3333333355553333333333111511e44433333333355533c3333333355333e444eeeeeeeeee44455f44f8ff444444e444eeeeeeeeee44455f44f8ff4444444444
33333333333355333333331111114444333333333335553333333333553e144444444444e444445fee4ff44444ee4444444444444444445ff44ff44444444444
3333333333888883331311111115544e3333cc3333333353333333311151144e444444ee444444ffffff5f444e4444444444444444444444ffff5f4444444444
3333333338888883331111111144455455333cccc333888883333331111115e444444e44444888844fff55fee444444e44444444444888844fff55f444444444
3333333338888888553111122224ee4535533333333888888531311111114554444ee444444888884ffff5ff44444ee444444444444888884ffff5ff44444444
555333333999999335ee4422222ee444333533333338888888531111111ee4554ee444444e88888844fff5ff4444e444444444444488888844fff5ff44444444
335555333999999ee55544222224444433355333333999999e55e11122224444e4444444e449999944fff5fff4ee4444444444444449999944fff5fff4444444
333335533ffffff4444555428884444433333553333999999445544222224444444444ee444999994efff5f4fe444444444444444449999944fff5f4f4444444
333eeee5eefffff444424458888844ee333eeee55eeffffff4444552222244ee44444e44444fffffe4ff55fef444444444444444444fffff44ff55f4f4444444
3eee444f54bbbb4444e22428888844e33eee44444544fffff4eee454288844e3444ee444444ffff444115114f444444444444444444ffff444115114f4444444
e444444f45bbbbbeeee42228888844e3e44444444554bbbb4ee44445888884e34ee4444444bbbbb444af5aa4f444444e4444444444bbbbb444af5aa4f4444444
4444eeeffbbbbbd444444428888844e34444eeeee445bbbbb4444442888884e3e44444444bbbbbbf4ffffaa444444ee4444444444bbbbbbf4ffffaa444444444
4eeee44ffdddddd444444448888844e34eeee444444bbbbbd4444422888882e34444444ffbbbbbbfffaaaaa44444e4444444444ffbbbbbbfffaaaaa444444444
ee44444ffbbbbbb444444448811144e3ee444444444ddddddf444442888882e3444444fffbbbbbb4eaa44aa444ee4444444444ff0bbbbbb44aa00aa444444444
444444444bbbbdd4444444411ddd44e344444444444bbbbbbff44444888884e344444ff4bbbbbbee4aa44aa44e44444444444ff4bbbbbb444aa00aa444444444
444444444dddddb44444444ddddd44e344444444444bbbbdd45ff444881114e34444ff44bbbbbb444aa44aaee44444444444ff44bbbbbb444aa00aa444444444
4444444441111114444444dddddd44ee44444444444dddddbfff444411ddd4ee444e4f44bbbbbb444aa44aa44444444444444f44bbbbbb444aa00aa444444444
44eee44448888884444444ddd422444444eee444444111111f444444ddddd4444ee44444111111444aa44aa44444444444444444111111444aa00aa444444444
eee3e444488888844444442244224444eee3e444444888888444444dddddd444e4444444888888444faeefa44444444e44444444888888444fa00fa444444444
e333e44448888ff44444442244224444e333e444444888888444444ddd4224444444444e88888844fff4fff444444ee44444444488888844fff0fff444444444
3333e4444ff44ff444444422444224443333e44444ff888ff44444422442244444444ee48888884fff4fff444444e444444444448888884fff4fff4444444444
33333e444ff44ff4444444222442224433333e444fff444ff4444442244224444444e444ff44ff4e4444444444ee444444444444ff00ff444444444444444444
33333e444ff444fff44444422244222433333e444ff4444ff44444422444224444ee444ff44ffee444444444ee4444444444444ff00fff444444444444444444
33333e444ff4444ff44444442244422433333e444fff4444fff44442224422244e4444fff4ffe4444444444e44444444444444fff0fff4444444444444444444
4444444444444444444444444444444444444a44644444dd44744448884444740000000000000000000000000000000000000000000000000000000000000000
ccccccccccccccccccccccccccccccccaccccaccc664444d44448888888847440000000000000000000000000000000000000000000000000000000000000000
44444444ddd444444444dddddddd4444aaa444aaaaa6444474488888888884480000000000000000000000000000000000000000000000000000000000000000
44444ddd44ddd444444dd444444ddd4444aa44aaaaaa644444888888888888880000000000000000000000000000000000000000000000000000000000000000
cccccdddd444ddccccc6666644444ddcccc4aaaaaaaaa64448888888888888840000000000000000000000000000000000000000000000000000000000000000
444ddd44ddd444d4444aaaa6666664d4444aaaaaaaaaa64448888888888888870000000000000000000000000000000000000000000000000000000000000000
44dd444444dd44d44aaaaaaaaaaa64d4aaaaaaaaaaaaa6444ffffffffffffff40000000000000000000000000000000000000000000000000000000000000000
44d44444444d4dd44aaaaaaaaaaaa6d44aaaaaaaaaaaa64d4ff777ff777ffff40000000000000000000000000000000000000000000000000000000000000000
ccddddddddddd4dccc666666666664dccc666666666664dc4ff7c7ff7c7ffff40000000000000000000000000000000000000000000000000000000000000000
44d44444444d44d444d44444444d44d444d44444444d44d44ff777ff777ffff40000000000000000000000000000000000000000000000000000000000000000
44d444dd444d44d444d444dd444d44d444d444dd444d44d44ffffffffffffff40000000000000000000000000000000000000000000000000000000000000000
44d44444444d44d444d44444444d44d444d44444444d44d44ffff7877fffff440000000000000000000000000000000000000000000000000000000000000000
ccd44444444d44dcccd44444444d44dcccd44444444d44dc44fff8888fffff440000000000000000000000000000000000000000000000000000000000000000
44d44444444d44d444d44444444d44d444d44444444d44d4444ff8888ffff4440000000000000000000000000000000000000000000000000000000000000000
44ddddddddddddd444ddddddddddddd444ddddddddddddd4444447778ffff4440000000000000000000000000000000000000000000000000000000000000000
444444444444444444444444444444444444444444444444444444fffffff4440000000000000000000000000000000000000000000000000000000000000000
__map__
0103070008000900000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

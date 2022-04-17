pico-8 cartridge // http://www.pico-8.com
version 35
__lua__

-- main
dt=1/30
g=-9.8/30 -- gravity for 1 frame
p={}
wind={}
waves={}
windtrails={}
cannonballs={}
next_windrail=0

speed_mult=0.1
max_speed=2
rot_speed=0.01
max_cannonball_dist=40
max_cannonball_halfdist=max_cannonball_dist * 0.5
max_wind=5
wave_sprites={4,5,6,5}
wave_sprites_count=count(wave_sprites)

function _init()
   p=make_boat(64, 64)
   wind={dir=rnd(1), str=rnd(max_wind)}
   for i=1,30 do
      add(waves, make_wave(64,64))
   end
end

function _update()
   update_wind()
   update_inputs()
   update_boat(p)
   update_cannonballs()
   update_waves()
   update_windtrails()
end

function _draw()
   cls(12)

   camera(p.x-64, p.y-64)
   draw_waves()
   draw_boat(p)
   draw_cannonballs()
   draw_windtrails()

   camera()
   print("wind: "..wind.dir.." "..wind.str)
   print("p: "..p.speed.." "..p.dir.." "..p.v)
   if count(cannonballs) > 0 then
      local ball = cannonballs[1]
      print("[ball] x: "..ball.x.." y: "..ball.y.." d: "..ball.dist.." state: "..ball.state)
   end
end

-- inputs and controls
function make_boat(x,y)
   return {x=x, y=y, speed=0, v=0, dir=0, rx=0, ry=0}
end
function update_inputs()
   if btn(0) then -- left
      p.dir += rot_speed
   elseif btn(1) then -- right
      p.dir -= rot_speed
   elseif btnp(2) then -- up
      p.speed += 1
   elseif btnp(3) then -- down
      p.speed -= 1
   end
   if btnp(4) then -- square
      make_cannonball(p, -1)
   end
   if btnp(5) then -- X
      make_cannonball(p, 1)
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
   p.rx, p.ry = dx - fdx, dy - fdy
   boat.x += fdx
   boat.y += fdy
end

boat_sprite_map = {x=0, y=0}
sail_sprites_map = {{x=2,y=0},{x=4,y=0},{x=6,y=0}}
function draw_boat(boat)
   draw_rotated_tile(boat.x, boat.y, -boat.dir, boat_sprite_map.x, boat_sprite_map.y, 1, false, 1.5)
   local sail_sprite_map = sail_sprites_map[boat.speed+1]
   draw_rotated_tile(boat.x, boat.y, -wind.dir, sail_sprite_map.x, sail_sprite_map.y, 1, false, 1.5)
end

-- states: [0: inair] [1: splash] [2: hit]
function make_cannonball(boat, side)
   local dir = normalize_angle(boat.dir+side*0.25)
   add(cannonballs, {x=boat.x, y=boat.y, dir=dir, v=1.0, dist=0, state=0, state_t=0})
end

function update_cannonballs()
   local toremove={}
   for i=1,count(cannonballs) do
      local ball = cannonballs[i]
      if ball.state == 0 then -- inair
         local dx, dy = ball.v*cos(ball.dir), ball.v*sin(ball.dir)
         -- todo: detect collision if z is withing detection threshold
         ball.x += dx
         ball.y += dy
         ball.dist += sqrt(dx*dx + dy*dy)
         if ball.dist > max_cannonball_dist then
            ball.state = 1
         end
      elseif ball.state == 1 then -- splash
         if ball.state_t < 30 then
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
         local r = lerp(2,5,invlerp(middist,max_cannonball_halfdist,0))
         circfill(ball.x, ball.y, r, 5)
      elseif ball.state == 1 then -- splash
         circfill(ball.x, ball.y, 2, 6) -- temp, to replace with sprite anim
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


__gfx__
000000000000000000000000000000000000000000000000000000000000000000000000000f0000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000770000000000000000000000000000000000000007f000000000000000000000000000000000000000000000000000
00700700445454400000000000000000070070000077000000770000000000000007000000077000000000000000000000000000000000000000000000000000
00077000444444440000000000000000000000000700700007007000000700000004700000047000000000000000000000000000000000000000000000000000
00077000444444440000000000000000000007700000077000000000000700000004700000047000000000000000000000000000000000000000000000000000
00700700445454400000000000000000000070070000700700000770000000000007000000077000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000700700000000000000000007f000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000f0000000000000000000000000000000000000000000000000000
__map__
0103070008000900000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

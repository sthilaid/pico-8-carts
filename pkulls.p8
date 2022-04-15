pico-8 cartridge // http://www.pico-8.com
version 35
__lua__

-- main
p={}
wind={}
waves={}
windtrails={}
next_windrail=0

speed_mult=0.1
max_speed=2
rot_speed=0.01
max_wind=5
wave_sprites={4,5,6,5}
wave_sprites_count=count(wave_sprites)

function _init()
   p=make_boat(64, 64)
   wind={dir=rnd(1), str=rnd(max_wind)}
   for i=1,30 do
      add(waves, make_wave())
   end
end

function _update()
   update_inputs()
   update_boat(p)
   update_wind()
   update_waves()
   update_windtrails()
end

function _draw()
   cls(12)
   draw_waves()
   draw_windtrails()
   draw_rotated_tile(p.x, p.y, -p.dir, 0, 0, 1, false, 1)
   draw_rotated_tile(p.x, p.y, -wind.dir, 2, 0, 1, false, 1)
   print("wind: "..wind.dir.." "..wind.str)
   print("p: "..p.speed.." "..p.dir.." "..p.v)
end

-->8
-- inputs and controls
function make_boat(x,y)
   return {x=x, y=y, speed=0, v=0, dir=0}
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
   p.dir = normalize_angle(p.dir)
   p.speed = min(max_speed, max(0, p.speed))
end

function update_boat(boat)
   delta_angle = shortest_angle_between_normalized(boat.dir, wind.dir)
   if delta_angle > 0.25 then
      boat.v = 0
   else
      boat.v = cos(delta_angle) * wind.str * boat.speed * speed_mult
   end
   boat.x += cos(boat.dir) * boat.v
   boat.y += sin(boat.dir) * boat.v
end

-->8
-- water and wind effects

function make_wave()
   speed = rnd(90)+60
   return {x=rnd(128), y=rnd(128), state=0, substate=rnd(speed), speed=speed, lifetime=rnd(15)}
end

function update_waves()
   for i=1,count(waves) do
      w = waves[i]
      w.substate += 1
      if w.substate >= w.speed then
         w.substate = 0
         w.state += 1
         if w.state >= w.lifetime then
            waves[i] = make_wave()
         end
      end
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
   add(windtrails, {x0=64+rnd(48)-24, y0=64+rnd(48)-24, life=0, length=0, maxlen=rnd(10)+3, maxdist=rnd(30)+15, dir=wind.dir})
end

function update_windtrails()
   next_windrail -= 1
   if next_windrail <=0 then
      make_windtrail()
      next_windrail = rnd(30)+30
   end

   toremove={}
   trailcount = count(windtrails)
   for i=1,trailcount do
      w = windtrails[i]
      if w.length >= w.maxdist then
         add(toremove, i)
      else
         w.life +=1
         mod = (max_wind - ceil(wind.str)) + 1
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
      l0 = max(0, w.length - w.maxlen)
      x0,y0 = w.x0 + cos(w.dir) * l0, w.y0 + sin(w.dir) * l0
      x1,y1 = w.x0 + cos(w.dir) * w.length, w.y0 + sin(w.dir) * w.length
      line(x0, y0, x1, y1, 6)
   end
end

-->8
-- utils
function normalize_angle(a)
   extra = flr(abs(a))
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
   delta = a2 - a1
   if delta < 0.5 then
      return delta
   else
      return 1.0 - delta
   end
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
0000000000000000000f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000007f00000000000007700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700445454400007700000000000070070000077000000770000000000000000000000000000000000000000000000000000000000000000000000000000
00077000444444440004700000000000000000000700700007007000000000000000000000000000000000000000000000000000000000000000000000000000
00077000444444440004700000000000000007700000077000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700445454400007700000000000000070070000700700000770000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000007f00000000000000000000000000000007007000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__map__
0103020000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

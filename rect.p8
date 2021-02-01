pico-8 cartridge // http://www.pico-8.com
version 29
__lua__

function rect(x, y, hw, hh, vx, vy, vangle, col)
   return {x=x, y=y, theta=0, hw=hw, hh=hh, vx=vx, vy=vy, vangle=vangle, col=col}
end

function updateRect(r)
   if ((r.x + r.vx > 128) or (r.x + r.vx < 0)) then
      r.vx = -r.vx
   end
   if ((r.y + r.vy > 128) or (r.y + r.vy < 0)) then
      r.vy = -r.vy
   end
   r.x += r.vx
   r.y += r.vy
   r.theta += r.vangle
   if r.theta > 1.0 then
      r.theta -= 1.0
   end
end

function drawRect(r)
   p0 = add(r, rotate(-r.hw, r.hh, r.theta))
   p1 = add(r, rotate(r.hw, r.hh, r.theta))
   p2 = add(r, rotate(r.hw, -r.hh, r.theta))
   p3 = add(r, rotate(-r.hw, -r.hh, r.theta))
   line(p0.x, p0.y, p1.x, p1.y, r.col)
   line(p2.x, p2.y, r.col)
   line(p3.x, p3.y, r.col)
   line(p0.x, p0.y, r.col)
end

function point(x, y)
   return {x=x, y=y}
end

function rotate(x, y, angle)
   theta = atan2(x, y) + angle
   norm = sqrt(x*x+y*y)
   return {x=cos(theta) * norm,
           y=sin(theta) * norm}
end

function add(v1, v2)
   return {x=v1.x + v2.x,
           y=v1.y + v2.y}
end

rects = {}

function _init()
   for i=1,100 do
      rects[i] = rect(rnd(128), rnd(128), 1+rnd(10), 1+rnd(10), rnd(5), rnd(5), rnd(0.2), flr(rnd(16)))
   end
end

function _update()
   for i=1,#rects do
      updateRect(rects[i])
   end
end

function _draw()
   cls()
   for i=1,#rects do
      drawRect(rects[i])
   end
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

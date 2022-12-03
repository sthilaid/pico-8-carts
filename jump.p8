pico-8 cartridge // http://www.pico-8.com
version 38
__lua__

p={}
ground=119
debug={}

function _init()
   p.x = 10
   p.y = 10
   p.vx = 0
   p.vy = 0
   p.w = 3
   p.h = 7
   p.state = 0 -- 0: ground, 1: inair
end

function _update()
   updatePlayerInputs()
   movePlayer()
end

function updatePlayerInputs()
   local dv = 1
   if (btn(0)) p.vx -= dv
   if (btn(1)) p.vx += dv
   if (btn(4)) updateJumpInput()
end

function updateJumpInput()
   if p.state == 0 and btnp(4) then
      p.vy -= 8
      p.state = 1
   end
end

function intersect(x,y,w,h,flag)
   for i=x,x+w do
      local cx = i \ 8
      for j=y,y+h do
         local cy = j \ 8
         local tile = mget(cx, cy)
         if fget(tile, flag) then
            return true
         end
      end
   end
   return false
end

function move(dx,dy)
   local safeX, safeY = p.x, safeY
   local dirX, dirY = dx / abs(dx), dy / abs(dy)
   for deltaX=0,abs(dx) do
      for deltaY=0,abs(dy) do
         local x, y = p.x + dirX*deltaX, p.y + dirY*deltaY
         local isColliding = intersect(x, y, p.w, p.h, 0)
         if isColliding then
            --stop()
            p.x, p.y = safeX, safeY
            return false
         else
            safeX, safeY = x, y
         end
      end
   end
   p.x = safeX
   p.y = safeY
   return true
end

function detectGround()
   return intersect(p.x, p.y+p.h+1, p.w, 0, 0)
end

function movePlayer()
   move(p.vx, p.vy)
   local onGround = detectGround()
   if onGround then
      p.vy = 0
      p.state = 0
   else
      p.vy += 1
      p.state = 1
   end
   p.vx = 0
end

function _draw()
   cls(0)
   camera(p.x-32,0)
   drawWorld()
   drawPlayer()
   camera()
   print("p: {"..p.x..","..p.y.."} v: {"..p.vx..","..p.vy.."}")
end

function drawPlayer()
   local x0 = p.x
   local x1 = x0 + p.w
   local y0 = p.y
   local y1 = y0 + p.h
   local c = 11
   if (p.state == 1) c = 14
   rectfill(x0,y0,x1,y1,c)
end

function drawWorld()
   map(0,0,0,0,128,64)
end

__gfx__
00000000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000cccccccc0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__gff__
0001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__map__
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000010100000000000000000000000000000000000000000000000000000001010000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000001010000000001010000000000000000000001010000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000010100000000000000000001010000000000000000000001010000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000010100000000000000010100000000000000000001010000000000000000000001010000000000000000000001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

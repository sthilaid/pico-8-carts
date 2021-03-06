pico-8 cartridge // http://www.pico-8.com
version 30
__lua__
--golf by sthilaid

-- global constants
btn_left    = 0
btn_right   = 1
btn_up      = 2
btn_down    = 3
btn_o       = 4
btn_x       = 5

flag_fairway= 0
flag_sand   = 1
flag_turf   = 2
flag_start  = 6
flag_hole   = 7

flag_green_right    = 1
flag_green_up       = 2
flag_green_left     = 3
flag_green_down     = 4

state_idle              = 0
state_swing_power       = 1
state_swing_accuracy    = 2
state_airborn           = 3
state_oob               = 4
state_green             = 5
state_rolling           = 6

g_gravity               = -9.8
g_aim_rot_accel         = 0.1
g_power_accel           = 0.5
g_accuracy_accel        = 0.5
g_zoom_accel            = 0.3
g_swing_power_target    = 0.85
g_swing_accuracy_target = 0.15
g_min_zoom_ratio        = 0.1
g_max_zoom_ratio        = 3

-- stats from https://blog.trackmangolf.com/trackman-average-tour-stats/
function make_club(name, maxdist, maxheight, powerTarget, powerTargetRange, accTarget, accTargetRange)
   g_swing_power_target    = 0.85
   g_swing_accuracy_target = 0.15
   return {name=name,maxdist=maxdist,maxheight=maxheight,
           powerTarget=powerTarget, powerTargetRange=powerTargetRange, accTarget=accTarget, accTargetRange=accTargetRange}
end
g_clubs     = { make_club("driver", 275, 32, 0.9, 0.01, 0.1, 0.01),
                make_club("3 wood", 243, 30, 0.85, 0.02, 0.15, 0.01),
                make_club("5 wood", 230, 31, 0.95, 0.02, 0.15, 0.02),
                make_club("3 iron", 212, 27, 0.8, 0.05, 0.2, 0.03),
                make_club("4 iron", 203, 28, 0.8, 0.05, 0.2, 0.03),
                make_club("5 iron", 194, 31, 0.8, 0.05, 0.2, 0.03),
                make_club("6 iron", 183, 30, 0.8, 0.05, 0.2, 0.03),
                make_club("7 iron", 172, 32, 0.8, 0.05, 0.2, 0.03),
                make_club("8 iron", 160, 31, 0.8, 0.05, 0.2, 0.03),
                make_club("9 iron", 148, 30, 0.8, 0.05, 0.2, 0.03),
                make_club("pw", 100, 29, 0.8, 0.1, 0.1, 0.04),
}
g_putter    = make_club("putter", 50, 0, 1, 0, 0, 1)

-- global state variables
g_lastFrameTime = 0
g_course        = false
g_state         = false
g_stateTime     = 0
g_swing_power   = 0
g_swing_accuracy= 0
g_aim_angle     = 0
g_club_index    = 1
g_zoom_ratio    = 1.0

g_ball_x, g_ball_y, g_ball_z    = 0,0,0
g_ball_vx, g_ball_vy, g_ball_vz = 0,0,0
g_ball_takeoff_x, g_ball_takeoff_y = 0,0
g_ball_ground_flags = 1 << flag_fairway

g_frame_offset_x = 0
g_frame_offset_y = 0

-- objets
function make_course(x0,y0,w,h,world_h)
   return {x0=x0,y0=y0,w=w,h=h,
           worldPixelRatio = h*8/world_h,
           pixelToWorldRatio = world_h/(h*8),
}
end
function set_course(course)
   g_course = course

   g_state         = state_idle
   g_stateTime     = 0
   g_swing_power   = 0
   g_swing_accuracy= 0
   g_aim_angle = -0.25

   -- init ball start by finding start coord in map
   g_ball_x, g_ball_y = 0,0
   for y=course.y0,course.y0+course.h do
      for x=course.x0,course.x0+course.w do
         local isStart = fget(mget(x,y), flag_start)
         if (isStart) then
            g_ball_x, g_ball_y = x*8+4 * course.pixelToWorldRatio, y*8+2 * course.pixelToWorldRatio
            goto done
         end
      end
   end
   ::done::
end
-->8
-- updates
function _init()
   pal(14,129)
   set_course(make_course(0,0,16,32,250,500))
   g_state  = state_idle
end

function _update()
   local t = time()
   local dt = t - g_lastFrameTime
   g_lastFrameTime = t

   update_state_transitions()
   update_state(dt)
   g_stateTime += dt
end

function _draw()
   cls()
   update_frame_offset()
   drawcourse(g_course)
   drawball()
   drawhud()
   drawdebug()
end

-->8
-- input & states

function calcSwingPower(clubdata)
   if (clubdata.powerTarget < 0) return g_swing_power

   local target, range = clubdata.powerTarget, clubdata.powerTargetRange
   local delta = g_swing_power - target
   if (abs(delta) < range) return 1.0
   local outrangeRatio = clamp(0,1,invlerp(0, 1-range, g_swing_power))
   return lerp(0.1, 1.0, outrangeRatio) -- should add min ratio (0.1) to club data?
end

function calcSwingAngle(clubdata)
   local target,range = clubdata.accTarget, clubdata.accTargetRange
   local delta = g_swing_accuracy - target
   if (abs(delta) < range) return g_aim_angle
   local total_range = g_swing_accuracy - (target + range)
   local outrangeRatio = clamp(0,1,invlerp(0, total_range, abs(delta) - range))
   local precisionAnglePenality = lerp(0.001, 0.125, outrangeRatio)
   local precisionDir = sgn(delta) -- (rnd() > 0.5 and 1) or -1
   local finalAngle = standardize_angle(g_aim_angle + precisionDir * precisionAnglePenality)
   -- print(g_swing_accuracy,0,40)
   -- print(outrangeRatio)
   -- print(precisionAnglePenality)
   -- print(finalAngle)
   -- stop()
   return finalAngle
end

function get_ball_flags()
   local px,py = g_ball_x * g_course.worldPixelRatio, g_ball_y * g_course.worldPixelRatio
   local mx,my = g_course.x0 + px \ 8, g_course.y0 + py \ 8
   local isOOB = (mx <= g_course.x0) or (mx >= g_course.x0+g_course.w) or (my <= g_course.y0) or (my >= g_course.y0+g_course.h)
   if (isOOB) return 0
   return fget(mget(mx,my))
end

function is_green(flags)
   if (flags & 1<<flag_fairway == 0) return false
   return (flags & (1<<flag_green_right | 1<<flag_green_up | 1<<flag_green_left | 1<<flag_green_down)) != 0
end

function transitionTo(newstate)
   local prevState = g_state
   g_state = newstate
   g_stateTime = 0

   -- on enter state
   if newstate == state_idle or newstate == state_green then
      g_swing_power, g_swing_accuracy = 0,0

      if prevState == state_oob then
         g_ball_x, g_ball_y = g_ball_takeoff_x, g_ball_takeoff_y
         g_aim_angle = -0.25
      end
      
   elseif newstate == state_swing_accuracy then
      g_swing_accuracy = g_swing_power

   elseif newstate == state_airborn then
      local maxdist     = g_clubs[g_club_index].maxdist
      local maxheight   = g_clubs[g_club_index].maxheight
      --local powerRatio  = lerp(0.5, 1.0, invlerp(-0.4, 0, g_swing_power - g_swing_power_target))
      local powerRatio  = calcSwingPower(g_clubs[g_club_index])
      local shotAngle   = calcSwingAngle(g_clubs[g_club_index])
      local dist,height = powerRatio * maxdist, powerRatio * maxheight
      local airTime     = sqrt(-8 * height / g_gravity)
      local speed       = dist / airTime
      g_ball_vx, g_ball_vy = speed * cos(shotAngle), speed * sin(shotAngle)
      g_ball_vz         = sqrt(-2*height*g_gravity)
      g_ball_takeoff_x, g_ball_takeoff_y = g_ball_x, g_ball_y
      -- print("pos: "..g_ball_x..","..g_ball_y..","..g_ball_z)
      -- print("vel: "..g_ball_vx..","..g_ball_vy..","..g_ball_vz)
      --print("powerRatio:"..powerRatio,0,4)
      -- print("d:".."h:"..height)
      -- print("airtime:"..airTime)
      -- print("vz: "..vz.." | "..g_ball_vz)
      --stop()
   elseif newstate == state_rolling then
      local maxspeed     = g_putter.maxdist / 5 -- dist over 5 secs
      local powerRatio  = calcSwingPower(g_putter)
      local shotAngle   = calcSwingAngle(g_putter)
      local speed        = powerRatio * maxspeed
      g_ball_vx, g_ball_vy = speed * cos(shotAngle), speed * sin(shotAngle)
      g_ball_vz         = 0
      g_ball_takeoff_x, g_ball_takeoff_y = g_ball_x, g_ball_y
      -- print("vx: "..g_ball_vx.." vy: "..g_ball_vy)
      -- stop()
   end
end

function update_state_transitions(dt)
   if g_state == state_idle or g_state == state_green then
      if (btnp(btn_x)) transitionTo(state_swing_power)
   elseif g_state == state_swing_power then
      if not btn(btn_x) then
         if is_green(g_ball_ground_flags) then  transitionTo(state_rolling)
         else                                   transitionTo(state_swing_accuracy)
         end
      end
   elseif g_state == state_swing_accuracy then
      if btn(btn_x) then
         transitionTo(state_airborn)
      end
   elseif g_state == state_oob then
      if (g_stateTime > 3) transitionTo(state_idle)
   end
end

function update_state(dt)
   if g_state == state_idle or g_state == state_green then
      local dir = 0
      if (btn(btn_right))   dir = 1
      if (btn(btn_left))    dir = -1
      if dir != 0 then
         g_aim_angle = standardize_angle(g_aim_angle + dir * g_aim_rot_accel * dt)
      end
      if btn(btn_o) then
         if (btn(btn_up))     g_zoom_ratio = clamp(g_min_zoom_ratio, g_max_zoom_ratio, g_zoom_ratio + g_zoom_accel * dt)
         if (btn(btn_down))   g_zoom_ratio = clamp(g_min_zoom_ratio, g_max_zoom_ratio, g_zoom_ratio - g_zoom_accel * dt)
      else
         if (btnp(btn_up))    g_club_index = clamp(1,#g_clubs, g_club_index-1)
         if (btnp(btn_down))  g_club_index = clamp(1,#g_clubs, g_club_index+1)
      end

   elseif g_state == state_swing_power then
      if btn(btn_x) then
         g_swing_power = clamp(0,1, g_swing_power + g_power_accel * dt)
      end

   elseif g_state == state_swing_accuracy then
      g_swing_accuracy = clamp(0,1, g_swing_accuracy - g_accuracy_accel * dt)

   elseif g_state == state_airborn then
      g_ball_x, g_ball_y = g_ball_x + g_ball_vx * dt, g_ball_y + g_ball_vy * dt
      g_ball_z = g_ball_z + g_ball_vz * dt

      if g_ball_z < 0 then
         g_ball_z = 0
         g_ball_vx, g_ball_vy, g_ball_vz = 0,0,0 -- no bounces for now
         g_ball_ground_flags = get_ball_flags()
         if g_ball_ground_flags == 0 then           transitionTo(state_oob)
         elseif is_green(g_ball_ground_flags ) then transitionTo(state_green)
         else                                       transitionTo(state_idle)
         end
      else
         g_ball_vz = g_ball_vz + g_gravity * dt
      end
   elseif g_state == state_rolling then
      g_ball_x, g_ball_y = g_ball_x + g_ball_vx * dt, g_ball_y + g_ball_vy * dt

      local sqdist = g_ball_vx * g_ball_vx + g_ball_vy * g_ball_vy
      if sqdist < 0.1 then
         g_ball_vx, g_ball_vy = 0,0
         g_ball_ground_flags = get_ball_flags()
         if g_ball_ground_flags == 0 then           transitionTo(state_oob)
         elseif is_green(g_ball_ground_flags ) then transitionTo(state_green)
         else                                       transitionTo(state_idle)
         end
      else
         local frictionDecay,minFrictionSpeed = 0.1, 0.1
         local sx,sy = sgn(g_ball_vx), sgn(g_ball_vy)
         g_ball_vx = g_ball_vx + (-1*sx)* max(minFrictionSpeed, frictionDecay * g_ball_vx)
         g_ball_vy = g_ball_vy + (-1*sy)* max(minFrictionSpeed, frictionDecay * g_ball_vy)
         if (sx != sgn(g_ball_vx)) g_ball_vx = 0
         if (sy != sgn(g_ball_vy)) g_ball_vy = 0
      end
   end
end

-->8
-- draw
function getZoom()
   local isOnGreen = is_green(g_ball_ground_flags)
   return (isOnGreen and g_max_zoom_ratio) or g_zoom_ratio
end
function worldToPixel(x,y)
   local zoomedPixelRatio = g_course.worldPixelRatio * getZoom()
   local wx,wy = x * zoomedPixelRatio - g_frame_offset_x, y * zoomedPixelRatio - g_frame_offset_y
   return wx,wy
end

function update_frame_offset()
   local zoom = getZoom()
   local zoomedPixelratio = g_course.worldPixelRatio * zoom
   local ballPX, ballPY = g_ball_x * zoomedPixelratio, g_ball_y * zoomedPixelratio
   local cx, cy = 64 + cos(g_aim_angle+0.5) * 42, 64 + sin(g_aim_angle+0.5) * 42
   g_frame_offset_x = clamp(0, max(0,g_course.w*8*zoom - 128), ballPX - cx)
   g_frame_offset_y = clamp(0, max(0,g_course.h*8*zoom - 128), ballPY - cy)
end

function drawcourse(course)
   local isGreen = g_state == state_green
   local zoom = getZoom()
   local pixelsPerSprite = 8 * zoom
   local offset_x, offset_y = g_frame_offset_x / pixelsPerSprite, g_frame_offset_y / pixelsPerSprite
   local width  = min(127, g_course.w * pixelsPerSprite - g_frame_offset_x)
   local height = min(127, g_course.h * pixelsPerSprite - g_frame_offset_y)
   for y=0,height do
      tline(0,y, width,y, offset_x, offset_y+y/pixelsPerSprite, 1/pixelsPerSprite, 0)
   end
   --print(offset_x.." | "..offset_y)
end

function drawball()
   local ballx, bally = worldToPixel(g_ball_x, g_ball_y)
   local ballr = lerp(1,8,invlerp(0,35, g_ball_z))
   circfill(ballx, bally, ballr, 7)
end

function drawClubHud()
   line(0,120,127,120,6)
   rectfill(0,121,127,127,14)
   local isgreen = is_green(g_ball_ground_flags)
   local club = (isgreen and g_putter) or g_clubs[g_club_index]
   local clubstr = "club: "..club.name
   color(5)
   -- print(clubstr,2,122)
   color(7)
   print(clubstr,1,122)
end
function drawAimHud()
   local pixeldist = g_clubs[g_club_index].maxdist * g_course.worldPixelRatio * getZoom()
   local ballx, bally = worldToPixel(g_ball_x, g_ball_y)
   local x0,y0 = ballx, bally
   local x1,y1 = x0 + cos(g_aim_angle) * pixeldist, y0 + sin(g_aim_angle) * pixeldist
   line(x0,y0, x1,y1, 7)
end
function drawhud()
   drawClubHud()

   if g_state == state_idle or g_state == state_green then
      drawAimHud()

   elseif g_state == state_swing_power or g_state == state_swing_accuracy then
      local x0,y0,w,h = 115,64,8,48
      local currentPad,targetPad = 3,1
      local powerHeight = h * (1.0 - g_swing_power)

      rect(x0-1,y0-1,x0+w+1,y0+h+1, 7) -- highlight
      rectfill(x0,y0,x0+w,y0+h, 6)
      
      if g_state == state_swing_power then
         local powerrange = g_clubs[g_club_index].powerTargetRange
         local targetRectHeight = h * (1.0-g_swing_power_target)
         local targetRectHeightRangeDelta = h * powerrange
         --line(x0-targetPad, y0+targetRectHeight, x0+w+targetPad, y0+targetRectHeight, 9)
         rectfill(x0, max(y0, y0+targetRectHeight-targetRectHeightRangeDelta),
                  x0+w, y0+targetRectHeight+targetRectHeightRangeDelta, 9)
         line(x0-currentPad, y0+powerHeight, x0+w+currentPad, y0+powerHeight,8)
         
         rectfill(x0,y0+powerHeight,x0+w,y0+h*g_swing_power+powerHeight,8)
      elseif g_state == state_swing_accuracy then
         local accuracyheight = h * (1 - g_swing_accuracy_target)
         local accuracyRange    = g_clubs[g_club_index].accTargetRange
         local targetRectHeight = h * (1.0-g_swing_accuracy_target)
         local targetRectHeightRangeDelta = h * accuracyRange
         line(x0-currentPad, y0+powerHeight, x0+w+currentPad, y0+powerHeight, 8)
         --line(x0-targetPad, y0+accuracyheight, x0+w+targetPad, y0+accuracyheight,9)
         rectfill(x0, max(y0, y0+targetRectHeight-targetRectHeightRangeDelta),
                  x0+w, y0+targetRectHeight+targetRectHeightRangeDelta, 9)
         
         local accuracyRatio = 1.0 - g_swing_accuracy
         line(x0-currentPad,y0+h*accuracyRatio,x0+w+currentPad,y0+h*accuracyRatio,10)
      end
   elseif g_state == state_oob then
      color(7)
      print("o.o.b.", 64,64)
      
   end
end

function drawdebug()
   cursor(0,0,7)
   print("s:"..g_state.." d:"..g_aim_angle.." p:"..g_swing_power.." a:"..g_swing_accuracy)
   print("vx:"..g_ball_vx.." vy:"..g_ball_vy)
   -- print("pos: "..format(g_ball_x)..","..format(g_ball_y)..","..format(g_ball_z))
   -- print("vel: "..format(g_ball_vx)..","..format(g_ball_vy)..","..format(g_ball_vz))
   print("flags: "..g_ball_ground_flags)
end

-->8
-- utils

function pow(base,exp)
   if exp <= 0 then return 1 end
   return base*pow(base,exp-1)
end

function format(num, digits)
   digits = digits or 2
   local v = pow(10, digits)
   return flr(num * v) / v
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

function clamp(minv,maxv,v)
   if (v < minv) return minv
   if (v > maxv) return maxv
   return v
end
function standardize_angle(angle)
   if (angle < -0.5) return angle+1
   if (angle > 0.5) return angle-1
   return angle
end
function lerp(a,b,ratio) return (b-a)*ratio + a end
function invlerp(a,b,v) return (v-a) / (b-a) end
function pbool(x) return x and "true" or "false" end

__gfx__
00000000bbbbbbbb333333333bbbbbbbbbbbbbb333333333aaaaaaaabbbbbb3339aaaaaaaaaaaa9333bbbbbb3333333339aaaaaaaaaaaaaaaaaaaa9333333333
00000000bbbbbbbb333333bb3bbbbbbbbbbbbbb3bb333333aaaaaaaabbbb339939aaaaaaaaaaa9939933bbbb9999999939aaaaaaaaaaaaaaaaaaaa9333333333
00700700bbbbbbbb3333bbbb33bbbbbbbbbbbb33bbbb3333aaaaaaaabbb399aab39aaaaaaaaaa93ba9993bbbaaaaaaaa39aaaaaaaaaaaaaaaaaaaa9333333333
00077000bbbbbbbb333bbbbb33bbbbbbbbbbbb33bbbbb333aaaaaaaabb399aaab399aaaaaaaa993baaa993bbaaaaaaaa39aaaaaaaaaaaaaaaaaaaa9333333333
00077000bbbbbbbb33bbbbbb333bbbbbbbbbb333bbbbbb33aaaaaaaab399aaaabb399aaaaaa993bbaaaa993baaaaaaaa39aaaaaaaaaaaaaaaaaaaa9333333333
00700700bbbbbbbb33bbbbbb3333bbbbbbbb3333bbbbbb33aaaaaaaab39aaaaabbb3999aaa993bbbaaaaa93baaaaaaaa39aaaaaaaaaaaaaaaaaaaa9333333333
00000000bbbbbbbb3bbbbbbb333333bbbb333333bbbbbbb3aaaaaaaa399aaaaabbbb33999933bbbbaaaaaa93aaaaaaaa39aaaaaa99999999aaaaaa9333333333
00000000bbbbbbbb3bbbbbbb3333333333333333bbbbbbb3aaaaaaaa39aaaaaabbbbbb3333bbbbbbaaaaaa93aaaaaaaa39aaaaaa33333333aaaaaa9333333333
bbbbbbbbbbbbbbbbbbbbbbbbb33333333333333bbbbbbbbb000000003333333339aaaaaaaaaaaa93333333333333333300000000000000000000000000000000
bbbbbbbbbbb88bbbbbbbbb33b33333333333333b33bbbbbb000000003333339939aaaaaaaaaaa993993333339999999900000000000000000000000000000000
bbbbbbbbbbb88bbbbbbb3333bb333333333333bb3333bbbb00000000333399aa339aaaaaaaaaa933a9993333aaaaaaaa00000000000000000000000000000000
bbbbbbbbbbb6bbbbbbb33333bb333333333333bb33333bbb0000000033399aaa3399aaaaaaaa9933aaa99333aaaaaaaa00000000000000000000000000000000
77bbbb77bbb6bbbbbb333333bbb3333333333bbb333333bb000000003399aaaa33399aaaaaa99333aaaa9933aaaaaaaa00000000000000000000000000000000
77bbbb77bb363bbbbb333333bbbb33333333bbbb333333bb00000000339aaaaa3333999aaa993333aaaaa933aaaaaaaa00000000000000000000000000000000
bbbbbbbbbb333bbbb3333333bbbbbb3333bbbbbb3333333b00000000399aaaaa3333339999333333aaaaaa93aaaaaaaa00000000000000000000000000000000
bbbbbbbbbbbbbbbbb3333333bbbbbbbbbbbbbbbb3333333b0000000039aaaaaa3333333333333333aaaaaa93aaaaaaaa00000000000000000000000000000000
00000000000000000000000000000000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
00000000000030000003003003000300bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
03003000000303000030030000303000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
00300300003000300300300000030000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
00030030000030000030030003000300bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
00300300000303000003003000303000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
03003000003000300000000000030000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb00000000000000000000000000000000000000000000000000000000
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
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f00000000000000000000000000000
__gff__
000101010101020202020202020202044181010101010002020202020000000000000000030509111f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__map__
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f02010101010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f01010101010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f01010110010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f01010101010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f03010101010101040f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f020101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f01010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f01010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0101010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0301010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f01010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0101010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
0f0f0f0301010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0101010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f170b1a0301010101010f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0c060e0f0101010101050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0c060e0f030101010101050f171a0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f180d190f0f0301010101010f0c061a0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f03010101010f0c06190f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f030101040f18190f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f022828050f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f011126260f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f032828040f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

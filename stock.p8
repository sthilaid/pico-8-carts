pico-8 cartridge // http://www.pico-8.com
version 38
__lua__

function makeCommodity(id,minp,maxp)
  return {id=id, minp=minp, maxp=maxp}
end

commodities = {
  makeCommodity("fe",400,800),
  makeCommodity("feo", 30, 100),
  makeCommodity("bbh", 1000, 3000)}

exchange={
  bids={},
  offers={}, 
  executed={}}

chartMode = 1
bookMode = 2

mode=bookMode
comIdx=1
barIdx=1
barOffset=0
period=40

function _init()
  for t=1,10000 do
    genorder(t)
  end
  menuitem(1, "Chart Mode", enableChartMode)
  menuitem(2, "Book Mode", enableBookMode)
end

function _update()
  updateInputs()
  updateTimeline()
end

function _draw()
  cls(0)

  local c = commodities[comIdx]
  if mode == chartMode then
    local lastT = exchange.executed[#exchange.executed].t
    drawCommodity(c.id, 0, lastT, period)
  elseif mode == bookMode then
    color(7)
    
    local orderedBids = findBook(c, exchange.bids, function (o1,o2) return o1.p > o2.p end)
    local orderedOffers = findBook(c, exchange.offers, function (o1,o2) return o1.p < o2.p end)
    local spread = 0
    if #orderedBids > 0 and #orderedOffers > 0 then
      spread = orderedOffers[1].p - orderedBids[1].p
    end
    print(c.id.." (spread: "..spread..")", 0, 0)
    color(9)
    displayBook("-- bids --", orderedBids, 0, 10)
    color(12)
    displayBook("-- offers --", orderedOffers, 64, 10)
  end
end

function findBook(c, book, comp)
  local orderedBook = {}
  for o in all(book) do
    if o.id == c.id then
      add(orderedBook, o)
    end
  end
  sort(orderedBook, comp)
  return orderedBook
end

function displayBook(title, book, xoffset, yoffset)
  local qtyoffset = 30
  print(title, xoffset, yoffset)
  yoffset += 10
  print("price", xoffset, yoffset)
  print("qty", xoffset+qtyoffset, yoffset)
  for o in all(book) do
    yoffset += 10
    print(o.p, xoffset, yoffset)
    print(o.qty, xoffset+qtyoffset, yoffset)
  end
end

function enableChartMode()
  mode = chartMode
end

function enableBookMode()
  mode = bookMode
end

function findmatch(order, tbl, op)
  for o in all(tbl) do
    if op(order, o) then
      return o
    end
  end
end

function order(id,p,qty,t,isbid)
  local bookTable;
  local matchTable;
  local op;

  if isbid then
    bookTable = exchange.bids
    matchTable = exchange.offers
    op = function (o1, o2) return o1.p >= o2.p end
  else
    bookTable = exchange.offers
    matchTable = exchange.bids
    op = function (o1, o2) return o1.p <= o2.p end
  end

  while qty > 0 do
    local neworder={id=id, p=p, qty=qty, t=t}
    local match = findmatch(neworder, matchTable, op)
    if match then
      neworder.p = match.p
      if match.qty > qty then
        match.qty = match.qty - qty
        add(exchange.executed, neworder)
        qty = 0
      else
        qty -= match.qty
        neworder.qty = match.qty
        add(exchange.executed, neworder)
        del(matchTable, match)
      end
    else
      add(bookTable, neworder)
      qty = 0
    end
  end
end

function genorder(t)
  local c = rnd(commodities)
  local price = c.minp + flr(rnd(c.maxp - c.minp))
  local qty = flr(rnd(1000))
  local isBid = rnd() < 0.5
  order(c.id, price, qty, t, isBid)
end

function updateInputs()
  if mode == chartMode then
    if btnp(2) then comIdx += 1 end
    if btnp(3) then comIdx -= 1 end
    if btnp(0) then barIdx -= 1 end
    if btnp(1) then barIdx += 1 end
    if btnp(5) then period *= 2 end
    if btnp(4) then period = max(10, flr(period/2)) end
  elseif mode == bookMode then
    if btnp(2) then comIdx += 1 end
    if btnp(3) then comIdx -= 1 end
  end
  if comIdx < 1 then comIdx = 1 end
  if comIdx > #commodities then comIdx = #commodities end
end

function updateTimeline()
end

bignum = 32000
function findPeriodData(id,from,to,period)
  local allPeriodData = {}
  local currentPeriod=from
  local nextPeriod=from+period
  local periodData={open=0, close=0, min=bignum, max=0, vol=0, count=0}
  for o in all(exchange.executed) do
    if o.t >= nextPeriod then
      add(allPeriodData, periodData)
      currentPrice = periodData.close
      periodData={open=currentPrice, close=currentPrice, min=currentPrice, max=currentPrice, vol=0, count=0}
      repeat 
        currentPeriod = nextPeriod
        nextPeriod = currentPeriod + period
        if o.t >= nextPeriod then add(allPeriodData, periodData) end
      until o.t < nextPeriod or currentPeriod > to

      periodData={open=currentPrice, close=currentPrice, min=bignum, max=0, vol=0, count=0}
    end
    if o.t > to then
      return allPeriodData
    end
    if o.t >= currentPeriod and o.t < nextPeriod and o.id == id then
      if periodData.open == 0 then periodData.open = o.p end
      periodData.close = o.p
      if o.p < periodData.min then periodData.min = o.p end
      if o.p > periodData.max then periodData.max = o.p end
      periodData.vol += o.qty
      periodData.count += 1
    end
  end
  return allPeriodData
end

function relerp(val, min, max, from, to)
  local ratio = (val-min) / (max-min)
  return from + (to-from)*ratio
end

function drawCommodity(id,from,to,period)
  local maxBid = 0
  for o in all(exchange.bids) do
    if o.p > maxBid then maxBid = o.p end
  end
  local minOffer = bignum
  for o in all(exchange.offers) do
    if o.p < minOffer then minOffer = o.p end
  end

  local pixelsPerBar = 3
  local maxBars = flr(128/pixelsPerBar)
  
  if barIdx < 1 then
    barOffset = max(0, barOffset - barIdx - 1)
    barIdx = 1
  end
  if barIdx > maxBars then
    barOffset = min(ceil(to/period)-maxBars, barOffset + (barIdx - maxBars))
    barIdx = maxBars
  end
  
  from = barOffset * period
  to = min(to, (barOffset + maxBars) * period)

  local periodDatas = findPeriodData(id,from,to,period)
  local totalMin = bignum
  local totalMax = 0
  for data in all(periodDatas) do
    if data.min < totalMin then totalMin = data.min end
    if data.max > totalMax then totalMax = data.max end
  end

  local barCount = min(maxBars, #periodDatas)
  for i=1,barCount do
    local data = periodDatas[i]
    local color = 7 --white
    if data.open < data.close then color = 11 end -- green
    if data.open > data.close then color = 8 end -- red
    if i==barIdx then color = 10 end -- yellow (selected)
    local bot = min(data.open, data.close)
    local top = max(data.open, data.close)
    local x0 = i*pixelsPerBar
    local y0 = relerp(bot, totalMin, totalMax, 100, 0)
    local x1 = x0+2
    local y1 = relerp(top, totalMin, totalMax, 100, 0)
    rectfill(x0, y0, x1, y1, color)

    if data.vol > 0 then
      local x0 = i*pixelsPerBar+1
      local y0 = relerp(data.min, totalMin, totalMax, 100, 0)
      local x1 = x0
      local y1 = relerp(data.max, totalMin, totalMax, 100, 0)
      rectfill(x0, y0, x1, y1, color)
    end
  end

  print("p: "..period, 0, 0, 7)
  print(from, 0, 94, 7)
  print(to, 110, 94, 7)
  rectfill(0, 100, 128, 128, 1)
  color(7)

  local selectedBar = periodDatas[barIdx]
  print(id.." vol: "..selectedBar.vol.."("..selectedBar.count..")", 1, 101)
  if selectedBar.vol > 0 then
    print("open: "..selectedBar.open.." close: "..selectedBar.close)
    print("min: "..selectedBar.min.." max: "..selectedBar.max)
  else
    print("price: "..selectedBar.close)
  end
  local selecctedTime = ((barOffset+barIdx-1))*period
  local selecctedTimeGlyphCount = 1
  if selecctedTime > 10 then selecctedTimeGlyphCount += 1 end
  if selecctedTime > 100 then selecctedTimeGlyphCount += 1 end
  if selecctedTime > 1000 then selecctedTimeGlyphCount += 1 end
  if selecctedTime > 10000 then selecctedTimeGlyphCount += 1 end
  print(selecctedTime, barIdx*pixelsPerBar+1-flr(selecctedTimeGlyphCount*4/2), 5, 7)
end

-------------------------------------------------------------------------------
-- utils

function sort(array, comp)
  function qsort(a, comp, lo, hi)
    function swap(a, i, j)
      local tmp = a[i];
      a[i] = a[j]
      a[j] = tmp
    end
    function part(a, comp, lo, hi)
      local pivot = a[hi]
      local i = lo-1
      for j = lo, hi-1 do
        if comp(a[j], pivot) then
          i += 1
          swap(a, i, j)
        end
      end
      i += 1
      swap(a, i, hi)
      return i
    end
    if lo >= hi or lo < 0 then return end
    local p = part(a, comp, lo, hi)
    qsort(a, comp, lo, p-1)
    qsort(a, comp, p+1, hi)
  end
  qsort(array, comp, 1, #array)
end

__gfx__
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00077000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00700700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

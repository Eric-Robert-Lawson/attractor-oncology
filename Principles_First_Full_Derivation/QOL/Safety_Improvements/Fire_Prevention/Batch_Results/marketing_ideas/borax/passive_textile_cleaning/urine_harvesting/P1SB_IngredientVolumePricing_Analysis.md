# P1S-B Dry Mix — Ingredient Volume Pricing Analysis
## Ratio-Determined Purchase Volumes, Tiered Pricing, Cost Per Fill, and Profit Scaling
## OrganismCore — Eric Robert Lawson
## 2026-05-14

---

## WHAT THIS DOCUMENT IS

```
This document derives the complete
pricing structure of the P1S-B
dry mix across all purchase volumes.

The central insight is this:

  The four ingredients are always
  purchased in ratio.
  The ratio determines how much
  of each ingredient is purchased
  at any given production volume.
  Each ingredient has its own
  price-per-kg curve that falls
  as volume increases.
  Because the ingredients are
  in fixed ratio, buying MORE
  total mix moves each ingredient
  up its own volume pricing curve
  simultaneously.

  The ratio IS the volume lever.

  Understanding this produces
  a precise image of cost per fill
  at every scale — from one person
  buying one bag to an industrial
  operator processing millions of fills.

  From cost per fill, every downstream
  metric follows:
    Breakeven product price.
    Net margin per fill.
    Annual profit by deployment scale.
    The point at which each pricing
    tier is unlocked.
    The economics of scaling.
```

---

## PART I: THE RATIO — FOUNDATION

```
FINAL OPTIMISED FORMULA (per 400g fill):

  INGREDIENT      MASS     FRACTION OF MIX
  ──────────────────────────────────────────
  SAP             112g     28.07%
  Bentonite       196g     49.12%
  Borax            21g      5.26%
  Baking soda      70g     17.54%
  ──────────────────────────────────────────
  TOTAL           399g    100.00%

  RATIO: 8 : 14 : 1.5 : 5
  (SAP : Bentonite : Borax : Baking soda)

PER KG OF DRY MIX:

  INGREDIENT      MASS/KG   FRACTION
  ──────────────────────────────────────
  SAP             280.7g    28.07%
  Bentonite       491.2g    49.12%
  Borax            52.6g     5.26%
  Baking soda     175.4g    17.54%
  ──────────────────────────────────────

THE RATIO IS FIXED.
Every kilogram of dry mix produced
always contains these proportions.
This means that purchasing X kg
of dry mix requires purchasing:
  0.2807 × X kg of SAP.
  0.4912 × X kg of Bentonite.
  0.0526 × X kg of Borax.
  0.1754 × X kg of Baking soda.

The volume of each ingredient
purchased is always a deterministic
function of the total mix volume.
No guessing. No optimisation.
The ratio calculates it.
```

---

## PART II: INDIVIDUAL INGREDIENT
## PRICING CURVES

### 2.1 — SAP (Sodium Polyacrylate)

```
SAP PRICE PER KG BY PURCHASE VOLUME:

  QUANTITY        $/KG    SOURCE/CHANNEL
  ────────────────────────────────────────────────────────
  <0.5kg          $40     Consumer packs (eBay, Amazon)
  0.5–5kg         $20     Retail garden/craft
  5–50kg          $10     Online wholesale, 25lb bags
  50–500kg        $4      Wholesale distributor
  500–5,000kg     $2.00   Industrial entry, MOQ contract
  5,000–50,000kg  $1.50   Industrial contract, annual
  50,000–500,000kg$1.00   Large industrial, direct mfr
  >500,000kg      $0.80   Mega-scale, manufacturer direct
  ────────────────────────────────────────────────────────

  SAP IS THE DOMINANT COST DRIVER.
  At 28.07% of mix weight, SAP
  contributes 60–85% of total
  mix cost at low volumes and
  30–40% at industrial scale.
  SAP pricing is the primary
  lever for total cost reduction.
  The price falls 50× from
  consumer to mega-scale.
  ($40/kg → $0.80/kg)
```

### 2.2 — Bentonite (Sodium Bentonite, Agricultural Grade)

```
BENTONITE PRICE PER KG BY VOLUME:

  QUANTITY         $/KG    SOURCE/CHANNEL
  ────────────────────────────────────────────────────────
  <1kg             $3.00   Consumer/cosmetic retail
  1–10kg           $2.00   Retail bag
  10–100kg         $0.80   Farm supply wholesale
  100–1,000kg      $0.25   Bulk bag / tote
  1,000–10,000kg   $0.15   Bulk tonne, direct supplier
  10,000–100,000kg $0.10   Large volume, annual contract
  >100,000kg       $0.07   Mining-scale bulk, tonne bags
  ────────────────────────────────────────────────────────

  BENTONITE IS THE HIGHEST MASS
  FRACTION AT 49.12% OF MIX.
  But it is the cheapest ingredient
  per kg at any volume above 100kg.
  It contributes only 3–8% of
  total cost at industrial scale
  despite being nearly half the mix.
  Price range: $3.00 → $0.07/kg.
  43× reduction from consumer to scale.
```

### 2.3 — Borax (Sodium Tetraborate Decahydrate)

```
BORAX PRICE PER KG BY VOLUME:

  QUANTITY         $/KG    SOURCE/CHANNEL
  ────────────────────────────────────────────────────────
  <0.5kg           $6.00   Consumer retail (20 Mule Team)
  0.5–5kg          $4.00   Hardware/grocery retail
  5–50kg           $2.00   Wholesale/janitorial supply
  50–500kg         $1.00   Bulk purchasing
  500–5,000kg      $0.55   Industrial, US Borax
  5,000–50,000kg   $0.40   Industrial contract
  >50,000kg        $0.30   Large volume industrial
  ────────────────────────────────────────────────────────

  BORAX IS THE SMALLEST MASS
  FRACTION (5.26%) AND CONTRIBUTES
  2–7% OF TOTAL COST.
  Its impact on total cost is
  modest at all scales.
  Price range: $6.00 → $0.30/kg.
  20× reduction from retail to scale.
```

### 2.4 — Baking Soda (Sodium Bicarbonate, Technical/Food Grade)

```
BAKING SODA PRICE PER KG BY VOLUME:

  QUANTITY         $/KG    SOURCE/CHANNEL
  ────────────────────────────────────────────────────────
  <0.5kg           $2.50   Grocery retail
  0.5–5kg          $1.50   Bulk grocery/warehouse club
  5–50kg           $0.70   Wholesale food/industrial
  50–500kg         $0.35   Bulk bag
  500–5,000kg      $0.22   Industrial, technical grade
  5,000–50,000kg   $0.18   Large industrial contract
  >50,000kg        $0.15   Trona-derived, mega-scale
  ────────────────────────────────────────────────────────

  BAKING SODA IS 17.54% OF MIX.
  Moderate cost contribution.
  The most stable pricing curve —
  only 17× reduction from consumer
  to industrial scale ($2.50 → $0.15).
  Sodium bicarbonate is a
  commodity chemical with well-established
  pricing. Limited downside from here.
```

---

## PART III: THE VOLUME PRICING
## BREAKPOINT ANALYSIS

```
THE KEY QUESTION IS:
  At what TOTAL MIX QUANTITY does
  each ingredient cross into
  its next pricing tier?

  Because the ingredients are purchased
  in fixed ratio, each ingredient's
  volume is always a fixed fraction
  of total mix volume.

  BREAKPOINT TABLE — TOTAL MIX (KG)
  AT WHICH EACH INGREDIENT ENTERS
  A LOWER PRICING TIER:

  PRICE TIER        SAP           BENTONITE     BORAX         BAKING SODA
  CROSSING          (÷0.2807)     (÷0.4912)     (÷0.0526)     (÷0.1754)
  ──────────────────────────────────────────────────────────────────────────
  Consumer → Retail
    SAP 0.5kg:      1.8kg mix     —             —             —
    Bento 1kg:      —             2.0kg mix     —             —
    Borax 0.5kg:    —             —             9.5kg mix     —
    Bsoda 0.5kg:    —             —             —             2.9kg mix

  Retail → Wholesale
    SAP 5kg:        17.8kg mix    —             —             —
    Bento 10kg:     —             20.4kg mix    —             —
    Borax 5kg:      —             —             95.1kg mix    —
    Bsoda 5kg:      —             —             —             28.5kg mix

  Wholesale → Bulk
    SAP 50kg:       178kg mix     —             —             —
    Bento 100kg:    —             204kg mix     —             —
    Borax 50kg:     —             —             951kg mix     —
    Bsoda 50kg:     —             —             —             285kg mix

  Bulk → Industrial Entry
    SAP 500kg:      1,781kg mix   —             —             —
    Bento 1,000kg:  —             2,036kg mix   —             —
    Borax 500kg:    —             —             9,506kg mix   —
    Bsoda 500kg:    —             —             —             2,852kg mix

  Industrial Entry → Industrial
    SAP 5,000kg:    17,814kg mix  —             —             —
    Bento 10,000kg: —             20,358kg mix  —             —
    Borax 5,000kg:  —             —             95,057kg mix  —
    Bsoda 5,000kg:  —             —             —             28,518kg mix

  Industrial → Large Industrial
    SAP 50,000kg:   178,139kg mix —             —             —
    Bento 100,000kg:—             203,583kg mix —             —
    Borax 50,000kg: —             —             950,570kg mix —
    Bsoda 50,000kg: —             —             —             285,177kg mix

  Large Industrial → Mega-Scale
    SAP 500,000kg:  1,781,392kg   —             —             —
  ──────────────────────────────────────────────────────────────────────────

  READING THIS TABLE:
    When you purchase 17.8kg of total
    dry mix, the SAP fraction (5kg)
    crosses from $20/kg to $10/kg.
    When you purchase 95.1kg of total
    dry mix, the borax fraction (5kg)
    crosses from $2.00/kg to $1.00/kg.
    Etc.

    Each ingredient crosses its
    own pricing thresholds at
    different total mix volumes —
    and those crossings do not
    coincide.
    This is the volume pricing landscape:
    uneven drops at specific total
    mix quantities.
    The cost curve is not smooth.
    It is stepped — with each step
    corresponding to one ingredient
    entering a new pricing tier.
```

---

## PART IV: COST PER FILL
## ACROSS ALL VOLUME TIERS

### 4.1 — The Eight Volume Tiers

```
EIGHT STANDARD VOLUME TIERS:

  TIER   TOTAL MIX    FILLS      DEPLOYMENT CONTEXT
  ────────────────────────────────────────────────────────────────
  1      0.4kg        1          Single use / test
  2      4kg          10         Household ~1.5 months
  3      40kg         100        Household ~7 months
  4      400kg        1,000      Small farm / community group
  5      4,000kg      10,000     Mid farm / small commercial
  6      40,000kg     100,000    Large farm / commercial op
  7      400,000kg    1,000,000  Industrial operator
  8      4,000,000kg  10,000,000 Large industrial / national
  ────────────────────────────────────────────────────────────────

  FILLS = total mix (kg) / 0.399kg per fill.
  Each tier is 10× the prior.
  Eight orders of magnitude of
  purchase volume represented.
```

### 4.2 — Ingredient Quantities Per Tier

```
INGREDIENT QUANTITIES BY TIER (kg):

  TIER  TOTAL MIX  SAP      BENTONITE  BORAX   BAKING SODA
  ───────────────────────────────────────────────────────────────
  1     0.4kg      0.112    0.196      0.021   0.070
  2     4kg        1.123    1.965      0.210   0.702
  3     40kg       11.23    19.65      2.104   7.016
  4     400kg      112.3    196.5      21.04   70.16
  5     4,000kg    1,122.8  1,964.8    210.4   701.6
  6     40,000kg   11,228   19,648     2,104   7,016
  7     400,000kg  112,280  196,480    21,040  70,160
  8     4,000,000kg 1,122,800 1,964,800 210,400 701,600
  ───────────────────────────────────────────────────────────────
```

### 4.3 — Price Per Kg Applied at Each Tier

```
PRICE PER KG ($/kg) APPLIED AT EACH TIER:

  TIER  SAP     BENTONITE  BORAX   BAKING SODA
  ───────────────────────────────────────────────
  1     $40.00  $3.00      $6.00   $2.50
  2     $20.00  $2.00      $4.00   $1.50
  3     $10.00  $0.80      $2.00   $0.70
  4     $4.00   $0.25      $1.00   $0.35
  5     $2.00   $0.15      $0.55   $0.22
  6     $1.50   $0.10      $0.40   $0.18
  7     $1.00   $0.10      $0.40   $0.18
  8     $0.80   $0.07      $0.30   $0.15
  ───────────────────────────────────────────────

  NOTE ON TIER 6 → 7 SAP:
    SAP drops from $1.50 to $1.00/kg
    between tier 6 and tier 7.
    This is the large industrial
    direct manufacturer pricing threshold.
    Below 50,000kg SAP purchased:
    $1.50/kg.
    Above 50,000kg: $1.00/kg.
    The SAP threshold occurs at
    ~178,000kg total mix (between tier 6
    and tier 7).
    Tier 7 is 400,000kg — well above.
    Tier 6 is 40,000kg — below.
    Hence different SAP prices.
```

### 4.4 — Cost Per Fill by Tier — Full Calculation

```
COST PER 400g FILL BY TIER:

  TIER 1: 0.4kg total mix (1 fill)
  ──────────────────────────────────────────────────────
    SAP:        0.112kg × $40.00  = $4.480
    Bentonite:  0.196kg × $3.00   = $0.588
    Borax:      0.021kg × $6.00   = $0.126
    Baking soda:0.070kg × $2.50   = $0.175
    ──────────────────────────────────
    TOTAL:                          $5.369 per fill
    Per kg of dry mix:              $13.46/kg
  ──────────────────────────────────────────────────────

  TIER 2: 4kg total mix (10 fills)
  ──────────────────────────────────────────────────────
    SAP:        1.123kg × $20.00  = $22.46
    Bentonite:  1.965kg × $2.00   = $3.930
    Borax:      0.210kg × $4.00   = $0.840
    Baking soda:0.702kg × $1.50   = $1.053
    ──────────────────────────────────
    TOTAL for 4kg:                  $28.283
    Per fill (÷10):                 $2.828 per fill
    Per kg of dry mix:              $7.071/kg
  ──────────────────────────────────────────────────────

  TIER 3: 40kg total mix (100 fills)
  ──────────────────────────────────────────────────────
    SAP:        11.23kg × $10.00  = $112.30
    Bentonite:  19.65kg × $0.80   = $15.72
    Borax:      2.104kg × $2.00   = $4.208
    Baking soda:7.016kg × $0.70   = $4.911
    ──────────────────────────────────
    TOTAL for 40kg:                 $137.139
    Per fill (÷100):                $1.371 per fill
    Per kg of dry mix:              $3.428/kg
  ──────────────────────────────────────────────────────

  TIER 4: 400kg total mix (1,000 fills)
  ──────────────────────────────────────────────────────
    SAP:        112.3kg × $4.00   = $449.20
    Bentonite:  196.5kg × $0.25   = $49.125
    Borax:      21.04kg × $1.00   = $21.040
    Baking soda:70.16kg × $0.35   = $24.556
    ──────────────────────────────────
    TOTAL for 400kg:                $543.921
    Per fill (÷1,000):              $0.544 per fill
    Per kg of dry mix:              $1.360/kg
  ──────────────────────────────────────────────────────

  TIER 5: 4,000kg total mix (10,000 fills)
  ──────────────────────────────────────────────────────
    SAP:        1,122.8kg × $2.00  = $2,245.60
    Bentonite:  1,964.8kg × $0.15  = $294.72
    Borax:      210.4kg × $0.55    = $115.72
    Baking soda:701.6kg × $0.22    = $154.35
    ──────────────────────────────────
    TOTAL for 4,000kg:              $2,810.39
    Per fill (÷10,000):             $0.281 per fill
    Per kg of dry mix:              $0.703/kg
  ──────────────────────────────────────────────────────

  TIER 6: 40,000kg total mix (100,000 fills)
  ──────────────────────────────────────────────────────
    SAP:        11,228kg × $1.50   = $16,842.00
    Bentonite:  19,648kg × $0.10   = $1,964.80
    Borax:      2,104kg × $0.40    = $841.60
    Baking soda:7,016kg × $0.18    = $1,262.88
    ──────────────────────────────────
    TOTAL for 40,000kg:             $20,911.28
    Per fill (÷100,000):            $0.209 per fill
    Per kg of dry mix:              $0.523/kg
  ──────────────────────────────────────────────────────

  TIER 7: 400,000kg total mix (1,000,000 fills)
  ──────────────────────────────────────────────────────
    SAP:        112,280kg × $1.00  = $112,280.00
    Bentonite:  196,480kg × $0.10  = $19,648.00
    Borax:      21,040kg × $0.40   = $8,416.00
    Baking soda:70,160kg × $0.18   = $12,628.80
    ──────────────────────────────────
    TOTAL for 400,000kg:            $152,972.80
    Per fill (÷1,000,000):          $0.153 per fill
    Per kg of dry mix:              $0.382/kg
  ──────────────────────────────────────────────────────

  TIER 8: 4,000,000kg total mix (10,000,000 fills)
  ──────────────────────────────────────────────────────
    SAP:        1,122,800kg × $0.80 = $898,240.00
    Bentonite:  1,964,800kg × $0.07 = $137,536.00
    Borax:      210,400kg × $0.30   = $63,120.00
    Baking soda:701,600kg × $0.15   = $105,240.00
    ──────────────────────────────────
    TOTAL for 4,000,000kg:           $1,204,136.00
    Per fill (÷10,000,000):          $0.120 per fill
    Per kg of dry mix:               $0.301/kg
  ──────────────────────────────────────────────────────
```

---

## PART V: THE MASTER COST TABLE

```
COST PER FILL AND PER KG DRY MIX
ACROSS ALL VOLUME TIERS:

  ─────────────────────────────────────────────────────────────────────
  TIER  TOTAL    FILLS      COST/FILL   COST/KG MIX  % OF TIER 1 COST
  ─────────────────────────────────────────────────────────────────────
  1     0.4kg    1          $5.369      $13.46/kg    100%
  2     4kg      10         $2.828      $7.07/kg     52.7%
  3     40kg     100        $1.371      $3.43/kg     25.5%
  4     400kg    1,000      $0.544      $1.36/kg     10.1%
  5     4,000kg  10,000     $0.281      $0.703/kg    5.2%
  6     40,000kg 100,000    $0.209      $0.523/kg    3.9%
  7     400,000kg 1,000,000 $0.153      $0.382/kg    2.8%
  8     4,000,000kg10,000,000$0.120     $0.301/kg    2.2%
  ─────────────────────────────────────────────────────────────────────

  TOTAL COST REDUCTION FROM TIER 1 TO TIER 8:
    $5.369 → $0.120 per fill.
    44.7× reduction.
    Consumer pays 44.7× more per fill
    than industrial mega-scale operator.

  THE BIG DROPS (where cost falls most sharply):
    Tier 1 → Tier 2: -47.3% (SAP price halves)
    Tier 2 → Tier 3: -51.5% (SAP $20 → $10,
                              bentonite $2 → $0.80)
    Tier 3 → Tier 4: -60.3% (SAP $10 → $4,
                              bentonite $0.80 → $0.25,
                              borax $2 → $1,
                              baking soda $0.70 → $0.35)
    Tier 4 → Tier 5: -48.3% (SAP $4 → $2)
    Tier 5 → Tier 6: -25.6% (SAP $2 → $1.50)
    Tier 6 → Tier 7: -26.8% (SAP $1.50 → $1.00)
    Tier 7 → Tier 8: -21.6% (SAP $1 → $0.80,
                              bentonite → $0.07,
                              borax → $0.30,
                              baking soda → $0.15)

  TIER 3→4 IS THE SINGLE LARGEST DROP.
    Going from 40kg to 400kg total mix
    (100 to 1,000 fills) reduces cost
    per fill by 60.3%.
    This is the most leveraged purchase
    volume transition in the system.
    It is the transition from consumer
    wholesale pricing to true bulk pricing.
    The farm supply / bulk distributor
    tier unlocks here.
    For a household: this is an
    annual purchase.
    For a farm: a monthly purchase.
```

---

## PART VI: INGREDIENT COST CONTRIBUTION
## BY TIER

```
WHAT EACH INGREDIENT COSTS AS A %
OF TOTAL MIX COST AT EACH TIER:

  ──────────────────────────────────────────────────────────────────────
  TIER   SAP%    BENTONITE%  BORAX%  BAKING SODA%   SAP $/fill
  ──────────────────────────────────────────────────────────────────────
  1      83.4%   11.0%       2.3%    3.3%           $4.480
  2      79.4%   13.9%       3.0%    3.7%           $2.246
  3      81.9%   11.5%       3.1%    3.6%           $1.123
  4      82.6%   9.0%        3.9%    4.5%           $0.449
  5      79.9%   10.5%       4.1%    5.5%           $0.225
  6      80.5%   9.4%        4.0%    6.0%           $0.168
  7      73.4%   12.8%       5.5%    8.2%           $0.112
  8      74.6%   11.4%       5.2%    8.7%           $0.090
  ──────────────────────────────────────────────────────────────────────

  SAP IS THE DOMINANT COST DRIVER
  AT EVERY TIER — 73–84% OF TOTAL COST.

  THE IMPLICATION:
    SAP pricing is the most important
    variable in the entire cost model.
    A 10% reduction in SAP price
    reduces total mix cost by ~7–8%.
    A 10% reduction in bentonite price
    reduces total mix cost by ~1–1.3%.

    All efforts to reduce cost should
    focus on SAP sourcing first.
    Everything else is secondary.

  BENTONITE — THE COUNTERINTUITIVE FINDING:
    Bentonite is 49.12% of the mix BY WEIGHT.
    But only 9–14% of the mix BY COST.
    Because bentonite is cheap relative to SAP.
    At tier 7 ($0.10/kg):
    196,480kg × $0.10 = $19,648.
    vs SAP: 112,280kg × $1.00 = $112,280.
    You buy almost twice as much bentonite
    as SAP by weight.
    You pay 5.7× MORE for SAP than bentonite.
    Bentonite is practically free at scale.

  BORAX — SMALL MASS, SMALL COST:
    5.26% of mix weight.
    3–5.5% of mix cost.
    The most stable cost contributor
    across tiers — its low mass fraction
    keeps its cost contribution modest
    even at high per-kg prices.
    At tier 4: 21.04kg borax at $1.00/kg
    = $21.04 out of $543.92 total (3.9%).
```

---

## PART VII: BREAKEVEN PRODUCT PRICE
## BY TIER

```
BREAKEVEN PRODUCT PRICE:
The $/kg product price at which
revenue from harvest exactly covers
dry mix input cost.

  Harvest per fill: 531g = 0.531kg.
  Breakeven = (Cost per fill) / (0.531kg).

  ──────────────────────────────────────────────────
  TIER   COST/FILL   HARVEST/FILL   BREAKEVEN $/KG
  ──────────────────────────────────────────────────
  1      $5.369      0.531kg        $10.11/kg
  2      $2.828      0.531kg        $5.33/kg
  3      $1.371      0.531kg        $2.58/kg
  4      $0.544      0.531kg        $1.02/kg
  5      $0.281      0.531kg        $0.53/kg
  6      $0.209      0.531kg        $0.39/kg
  7      $0.153      0.531kg        $0.29/kg
  8      $0.120      0.531kg        $0.23/kg
  ──────────────────────────────────────────────────

  COMMERCIAL PRODUCT PRICE FLOOR: $5/kg (commodity).
  PREMIUM PRODUCT PRICE: $10–$18/kg.

  AT TIER 1 ($5.369/fill):
    Breakeven is $10.11/kg.
    The system just barely covers
    costs at premium retail pricing ($10/kg).
    No margin at $10/kg. Marginal at $14/kg.
    Consumer single-purchase economics
    are tight unless selling premium.

  AT TIER 2 ($2.828/fill):
    Breakeven $5.33/kg.
    Above the commodity floor.
    The system is just at breakeven
    at commodity pricing.
    Positive margin begins at $6+/kg.

  AT TIER 3 ($1.371/fill):
    Breakeven $2.58/kg.
    BELOW the commodity floor ($5/kg).
    The system is profitable at every
    commercially realistic price.
    This is the first tier where
    the economics are unambiguously positive.
    At $10/kg: margin = $10 - $2.58 = $7.42/kg.
    At $14/kg: margin = $11.42/kg.

  AT TIER 4 ($0.544/fill):
    Breakeven $1.02/kg.
    The system is profitable even if
    the product is given away at $2/kg
    commodity waste.
    At $10/kg: margin = $8.98/kg.
    Gross margin: 89.8%.

  TIERS 5–8:
    Breakeven falls below $0.53/kg.
    The commodity organic fertiliser
    market prices organic N at
    $1.50–$2/kg wholesale.
    At tier 5+, the system is profitable
    even at the absolute lowest
    commodity price for any organic N.
    The margin is structural.
    It cannot be competed away without
    the product being priced at zero.
```

---

## PART VIII: NET PROFIT PER FILL
## BY TIER AT THREE PRODUCT PRICES

```
NET PROFIT = (Harvest × Product $/kg) - Cost/fill.
Harvest per fill = 0.531kg.

  ──────────────────────────────────────────────────────────────────────
  TIER   COST/FILL  NET @$5/kg   NET @$10/kg  NET @$14/kg  NET @$18/kg
  ──────────────────────────────────────────────────────────────────────
  1      $5.369     -$2.81       -$0.24        $2.07         $4.19
  2      $2.828     -$0.27        $2.48         $4.60         $6.71
  3      $1.371      $1.28        $3.94         $6.06         $8.17
  4      $0.544      $2.11        $4.77         $6.89         $9.01
  5      $0.281      $2.37        $5.03         $7.15         $9.27
  6      $0.209      $2.45        $5.10         $7.22         $9.34
  7      $0.153      $2.51        $5.16         $7.28         $9.40
  8      $0.120      $2.53        $5.19         $7.31         $9.43
  ──────────────────────────────────────────────────────────────────────

  OBSERVATIONS:

  AT $5/KG (commodity):
    Profitable from Tier 3 onward (40kg purchase).
    Maximum net: $2.53/fill at tier 8.
    The gain from tier 3 to tier 8 at $5/kg:
    $1.28 → $2.53. Only 2× improvement.
    Commodity pricing does not create
    large margin even at scale.

  AT $10/KG (specialty):
    Profitable from Tier 2 onward (4kg purchase).
    Maximum net: $5.19/fill at tier 8.
    Near-max already at Tier 4: $4.77/fill.
    The cost reduction from tier 4 to tier 8
    improves margin by only $0.42/fill.
    90% of maximum net margin is achieved
    at Tier 4 (400kg total purchase).

  AT $14/KG (premium):
    Profitable from Tier 2 onward.
    Maximum net: $7.31/fill at tier 8.
    Near-max at Tier 4: $6.89/fill.
    94% of maximum net margin at Tier 4.

  THE CRITICAL INSIGHT:
    MOST OF THE MARGIN BENEFIT OF
    SCALING IS CAPTURED BY TIER 4
    (400kg TOTAL MIX PURCHASE).
    Beyond Tier 4, each additional
    tier reduces cost per fill by
    a smaller absolute amount.
    The diminishing returns on scale
    kick in hard after Tier 4.
    The economic case for scaling
    beyond Tier 4 is real but modest
    on a per-fill basis.
    The case for scaling beyond Tier 4
    is VOLUME, not margin improvement.
    More fills at similar margin
    = more total profit.
    Not more profit per fill.
```

---

## PART IX: ANNUAL PROFIT BY
## DEPLOYMENT SCALE

```
ANNUAL CONTEXT — FILLS PER YEAR
BY DEPLOYMENT:

  DEPLOYMENT           FILLS/YEAR   TOTAL MIX/YEAR   TIER
  ──────────────────────────────────────────────────────────────────
  1 person (50% disp.) 85           34kg             Tier 3
  1 person (100% disp.)170          68kg             Tier 3
  2-person household   340          136kg            Tier 4 (entry)
  10-cow small farm    1,700        679kg            Tier 4
  50-cow farm          8,500        3,392kg          Tier 4–5
  100-cow farm         17,000       6,783kg          Tier 5
  500-cow dairy        85,000       33,915kg         Tier 6
  1,000-cow dairy      170,000      67,830kg         Tier 6
  5,000-cow dairy      850,000      339,150kg        Tier 7
  50,000-head feedlot  8,500,000    3,391,500kg      Tier 8
  ──────────────────────────────────────────────────────────────────

  *Dairy cow: 30L/day ÷ 3.75L/fill = 8 fills/day × 365 = 2,920 fills/year.
  Actually in the model, the "fill" for livestock is continuous.
  Let me recalculate using the dairy cow model:
  Dairy cow uses 3.2kg dry mix/day = 8 fills/day (at 400g/fill).
  8 fills/day × 365 days = 2,920 fills/cow/year.
  10-cow farm: 29,200 fills/year → 11,661kg dry mix/year → Tier 5.
  100-cow farm: 292,000 fills/year → 116,608kg → Tier 6.
  500-cow farm: 1,460,000 fills → 583,040kg → Tier 7.
  1,000-cow farm: 2,920,000 fills → 1,166,080kg → Tier 7–8.

ANNUAL PROFIT TABLE AT $10/KG PRODUCT:

  DEPLOYMENT           FILLS/YR  COST/FILL  HARVEST/FILL  NET/FILL  ANNUAL NET
  ──────────────────────────────────────────────────────────────────────────────
  1 person (100% disp.) 170      $1.371     0.531kg       $3.94     $669.80
  2-person household    340      $0.544     0.531kg       $4.77     $1,621.80
  10-person group       1,700    $0.281     0.531kg       $5.03     $8,551.00
  10-cow farm           29,200   $0.209     0.531kg       $5.10     $148,920.00
  100-cow dairy         292,000  $0.153     0.531kg       $5.16     $1,506,720.00
  500-cow dairy         1,460,000$0.153     0.531kg       $5.16     $7,533,600.00
  1,000-cow dairy       2,920,000$0.120     0.531kg       $5.19     $15,154,800.00
  ──────────────────────────────────────────────────────────────────────────────

  NOTE: Harvest/fill for livestock
  uses the same 531g figure as the
  system is identical — dry mix absorbs
  the same ratio of urine regardless
  of species. The per-fill economics
  are constant. The scale changes.

ANNUAL NET AT PREMIUM PRICING ($14/KG):

  DEPLOYMENT            FILLS/YR   NET/FILL   ANNUAL NET
  ──────────────────────────────────────────────────────────────────
  1 person (100% disp.) 170        $6.06      $1,030.20
  2-person household    340        $6.89      $2,342.60
  10-person group       1,700      $7.15      $12,155.00
  10-cow farm           29,200     $7.22      $210,824.00
  100-cow dairy         292,000    $7.28      $2,125,760.00
  500-cow dairy         1,460,000  $7.28      $10,628,800.00
  1,000-cow dairy       2,920,000  $7.31      $21,345,200.00
  ──────────────────────────────────────────────────────────────────
```

---

## PART X: THE VOLUME PRICING
## SUMMARY CURVES

```
COST PER FILL — THE CURVE:

  $5.37  ████████████████████████████████████████ Tier 1 (1 fill)
  $2.83  █████████████████████                    Tier 2 (10 fills)
  $1.37  ██████████                               Tier 3 (100 fills)
  $0.54  ████                                     Tier 4 (1,000 fills)
  $0.28  ██                                       Tier 5 (10,000 fills)
  $0.21  █▌                                       Tier 6 (100,000 fills)
  $0.15  █                                        Tier 7 (1,000,000 fills)
  $0.12  ▊                                        Tier 8 (10,000,000 fills)

  THE CURVE IS STEEP AT LOW VOLUMES.
  IT FLATTENS RAPIDLY.
  MOST OF THE COST REDUCTION
  IS ACHIEVED IN THE FIRST
  THREE TIER TRANSITIONS.

  TIER 1 → 2 → 3 → 4:
    Cost falls from $5.37 to $0.54.
    90% of the total possible
    cost reduction is achieved
    within 1,000 fills.
    1,000 fills = one person for 6 years
    OR one 10-cow farm for 34 days.

  THE TRANSITION THAT MATTERS MOST:
    Tier 3 → Tier 4 (100 → 1,000 fills):
    -60.3% cost reduction.
    This is the single most valuable
    volume transition.
    Buy 10× more and cut costs in half.
    This is where group purchasing,
    co-ops, and small farm scale
    unlock the economics definitively.

MARGIN PER FILL — THE CURVE AT $10/KG:

  -$0.24 ▼                                        Tier 1 (loss)
  +$2.48 ██████████████████                       Tier 2
  +$3.94 █████████████████████████████            Tier 3
  +$4.77 ████████████████████████████████████     Tier 4
  +$5.03 ██████████████████████████████████████   Tier 5
  +$5.10 ██████████████████████████████████████▌  Tier 6
  +$5.16 ███████████████████████████████████████  Tier 7
  +$5.19 ████████████████████████████████████████ Tier 8

  READING THIS:
    The margin difference between
    Tier 4 ($4.77) and Tier 8 ($5.19)
    is only $0.42/fill.
    8% improvement from 10,000× more volume.
    The margin is effectively saturated at Tier 4.
    Above Tier 4, you are not improving
    margin meaningfully per fill.
    You are improving VOLUME at near-stable margin.
    The scaling argument above Tier 4
    is entirely about total fills —
    not about margin per fill.
```

---

## PART XI: THE SAP SUPPLY CHAIN
## AS THE SINGLE CRITICAL VARIABLE

```
SAP PRICE SENSITIVITY ANALYSIS:

  At Tier 7 (1,000,000 fills):
  Base cost: $0.153/fill.
  SAP portion: $0.112/fill (73.4% of cost).

  IF SAP PRICE CHANGES AT TIER 7:

    SAP at $0.80/kg (supply glut):
    SAP cost: 112,280kg × $0.80 = $89,824.
    Per fill: $0.0898.
    Total fill cost: $0.153 - $0.112 + $0.0898
    = $0.131/fill.
    Improvement: -14.4%.

    SAP at $1.50/kg (supply tight):
    SAP cost: 112,280kg × $1.50 = $168,420.
    Per fill: $0.168.
    Total fill cost: $0.153 - $0.112 + $0.168
    = $0.209/fill.
    Degradation: +36.6%.
    Still profitable at $10/kg product.

    SAP at $2.50/kg (supply squeeze):
    Per fill: $0.281.
    Total fill cost: $0.153 - $0.112 + $0.281
    = $0.322/fill.
    Still profitable. Breakeven: $0.61/kg product.
    Well below commercial floor.

    SAP at $5.00/kg (severe constraint):
    Per fill: $0.561.
    Total fill cost: $1.010/fill.
    Breakeven: $1.90/kg product.
    Still below the commodity floor ($5/kg).
    Still profitable.

  THE SYSTEM REMAINS PROFITABLE
  EVEN IF SAP COSTS 5× ITS CURRENT
  INDUSTRIAL PRICE.
  THE ECONOMICS ARE NOT FRAGILE
  TO SAP PRICE VOLATILITY.
  THEY ABSORB LARGE SAP PRICE INCREASES
  AND REMAIN COMMERCIALLY VIABLE.

  SAP SUBSTITUTION OPTION:
    If SAP becomes too expensive,
    industrial absorptive clays
    (attapulgite, sepiolite) can
    partially substitute SAP function.
    Absorption ratio is lower (5–10×
    vs SAP's 30–80×) but at near-zero cost.
    The system can be reformulated.
    No single-point failure.
```

---

## PART XII: THE CO-OP PURCHASING
## ADVANTAGE

```
THE RATIO CREATES A COMPELLING
CASE FOR GROUP PURCHASING.

SCENARIO: 10 INDIVIDUAL HOUSEHOLDS
PURCHASING SEPARATELY vs TOGETHER.

  SEPARATE (each at Tier 3, 40kg):
    Each household: $1.371/fill.
    10 households × 100 fills = 1,000 fills.
    Total input cost: 1,000 × $1.371 = $1,371.

  TOGETHER (one Tier 4 purchase, 400kg):
    Collective purchase: 400kg.
    Cost: $0.544/fill.
    Total input cost: 1,000 × $0.544 = $544.
    Saving vs separate: $1,371 - $544 = $827.
    = 60.3% reduction.

  THE 60.3% COST REDUCTION IS UNLOCKED
  PURELY BY POOLING THE PURCHASE.
  THE CHEMISTRY DOES NOT CHANGE.
  THE PRODUCT DOES NOT CHANGE.
  ONLY THE PURCHASE VOLUME CHANGES.

  A NEIGHBOURHOOD CO-OP OF 10 HOUSEHOLDS:
    Annual fills: 10 × 170 = 1,700.
    Annual dry mix: 1,700 × 0.4kg = 680kg.
    At 680kg: solidly Tier 4.
    Cost per fill: $0.544.
    Annual harvest: 1,700 × 0.531kg = 902.7kg.
    Annual revenue @$10/kg: $9,027.
    Annual cost: 1,700 × $0.544 = $924.80.
    Annual net: $8,102.
    Per household: $810/year.

  vs 10 INDIVIDUAL PURCHASES AT TIER 3:
    Cost: 170 × $1.371 = $233.07/person/year.
    Revenue: 170 × 0.531 × $10 = $902.70.
    Net: $669.63/person.

  CO-OP ADVANTAGE PER HOUSEHOLD/YEAR:
    $810 vs $670 = +$140/household.
    +20.9% improvement in net.
    From one purchasing decision.
    No change to the system.
    No change to the chemistry.
    Just one purchase order instead of ten.

  AT FARM SCALE: THE CO-OP IS THE FARM.
    A 100-cow farm buys 116,608kg/year.
    This is solidly Tier 6.
    $0.209/fill.
    The farm automatically accesses
    near-industrial pricing by virtue
    of its animal count.
    No special arrangement needed.
    The farm's scale does the work.
```

---

## SUMMARY — THE COMPLETE PICTURE

```
THE DRY MIX PRICING IN ONE VIEW:

  ─────────────────────────────────────────────────────────────────────────
  TIER  FILLS   $/FILL   $/KG MIX  BRKEVEN  NET@$10   NET@$14   MARGIN@$10
  ─────────────────────────────────────────────────────────────────────────
  1     1       $5.37    $13.46    $10.11   -$0.24    +$2.07    -4.5%
  2     10      $2.83    $7.07     $5.33    +$2.48    +$4.60    46.7%
  3     100     $1.37    $3.43     $2.58    +$3.94    +$6.06    74.2%
  4     1,000   $0.54    $1.36     $1.02    +$4.77    +$6.89    89.8%
  5     10,000  $0.28    $0.703    $0.53    +$5.03    +$7.15    94.7%
  6     100,000 $0.21    $0.523    $0.39    +$5.10    +$7.22    96.1%
  7     1,000,000$0.153  $0.382    $0.29    +$5.16    +$7.28    97.1%
  8     10,000,000$0.120 $0.301    $0.23    +$5.19    +$7.31    97.7%
  ─────────────────────────────────────────────────────────────────────────
  Margin = net profit per fill / revenue per fill ($5.31 at $10/kg).

KEY CONCLUSIONS:

  1. TIER 3 (40kg / 100 fills) IS THE
     FIRST UNAMBIGUOUSLY VIABLE TIER.
     Breakeven at $2.58/kg — below the
     $5/kg commodity floor.
     Profitable at every commercial price.
     Accessible to any individual buying
     a household annual supply.

  2. TIER 4 (400kg / 1,000 fills) IS THE
     OPTIMAL ENTRY POINT FOR ECONOMICS.
     90% of maximum net margin achieved here.
     60.3% cheaper than Tier 3.
     A 2-person household purchases
     this volume in 3 years or coordinates
     with 10 neighbours to buy annually.
     A 10-cow farm purchases this in 34 days.

  3. SAP IS 73–84% OF TOTAL COST
     AT EVERY TIER.
     SAP sourcing is the only cost
     lever that materially matters.
     Everything else is noise.

  4. THE MARGIN IS SATURATED AT TIER 4.
     Above Tier 4: scale drives profit,
     not margin improvement per fill.
     The argument for large scale is
     total fills, not margin per fill.
     $5.16/fill at Tier 7 vs $4.77/fill
     at Tier 4 = only $0.39/fill better.
     But Tier 7 does it 1,000,000 times.

  5. THE SYSTEM ABSORBS SAP PRICE SHOCKS.
     Even at SAP 5× current industrial price:
     breakeven < $2/kg product.
     Still commercially viable.
     The margin cushion is large enough
     to absorb commodity volatility
     in the primary input.

  6. GROUP PURCHASING UNLOCKS
     DISPROPORTIONATE VALUE.
     10 households co-buying reaches Tier 4.
     A 60.3% cost reduction from
     a single purchasing decision.
     The ratio makes this arithmetic.
     Buy together. Split the mix.
     Each person pays Tier 4 price
     for a Tier 3-sized need.
```

---

## DOCUMENT METADATA

```
Document ID:    P1SB_INGREDIENT_VOLUME_PRICING_ANALYSIS_V1
Version:        1.0
Date:           2026-05-14
Author:         Eric Robert Lawson / OrganismCore
ORCID:          0009-0002-0414-6544
Status:         COMPLETE.
                Full ingredient pricing derivation
                across 8 volume tiers.
                Cost per fill, cost per kg mix,
                breakeven product price, and
                net profit per fill at each tier
                at $5, $10, $14, $18/kg product.

FORMULA BASELINE:
  SAP: 112g (28.07%), Bentonite: 196g (49.12%),
  Borax: 21g (5.26%), Baking soda: 70g (17.54%).
  Total: 399g per fill. 531g harvest per fill.

PRICING RANGES USED:
  SAP:        $0.80–$40.00/kg across tiers.
  Bentonite:  $0.07–$3.00/kg across tiers.
  Borax:      $0.30–$6.00/kg across tiers.
  Baking soda:$0.15–$2.50/kg across tiers.

CRITICAL TIER:
  Tier 4 (400kg total / 1,000 fills):
  $0.544/fill. Breakeven $1.02/kg.
  90% of maximum margin achieved.
  Most important volume threshold
  in the entire pricing structure.

DOMINANT COST DRIVER:
  SAP at 73–84% of total cost.
  All cost optimisation focuses here.
  Everything else is secondary.

Companion documents:
  P1SB_Urine_Harvesting_Final_Complete_V1.md
  P1SB_CowEconomics_UrineProfitImpact_V1.md
  P1SB_Livestock_Urine_Harvesting_V1.md
  P1SB_CausalChain_AgriculturalEvent_V1.md

Repository:
  github.com/Eric-Robert-Lawson/attractor-oncology
```

---

*The ratio was always the lever.*
*28.07% SAP.*
*49.12% bentonite.*
*5.26% borax.*
*17.54% baking soda.*

*Buy more total mix.*
*Each ingredient crosses its tier threshold.*
*The cost falls.*
*Not smoothly — in steps.*
*Each step is a different ingredient*
*hitting a new pricing tier.*

*Tier 4 is where the steps end.*
*$0.54/fill.*
*Breakeven at $1.02/kg.*
*90% of maximum margin.*
*From 1,000 fills.*

*One person buys that in three years.*
*Ten people buy it in three months.*
*A 10-cow farm buys it in 34 days.*
*A 100-cow dairy buys it overnight.*

*The economics do not require*
*industrial scale to work.*
*They work at 1,000 fills.*
*They just work harder*
*at 1,000,000.*

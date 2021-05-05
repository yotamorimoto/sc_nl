
Plogist : Pattern {
	var <>r, <>xi, <>length;

	*new { |r=1.5, xi=0.1, length=inf|
		^super.newCopyArgs(r, xi, length);
	}

	embedInStream { |inval|
		var xn = xi;
		var rStrm = r.asStream;
		var rVal;
		length.do({
			rVal = rStrm.next(inval);
			if(rVal.isNil, { ^inval });
			xn = 1.0 - (rVal * xn.squared);
			inval = xn.yield;
		});
		^inval;
	}
}

Pcml : ListPattern {
	var <>r, <>g;

	*new { |list, r=1.5, g=0.1, repeats=inf|
		^super.new(list, repeats).r_(r).g_(g);
	}

	f { |r,x| ^1.0 - (r * x.squared) }

	evolve {| prev, r, g |
		var next = Array.newClear(prev.size), halfG = g * 0.5;
		prev.size.do({|i|
			next[i] = ((1.0 - g) * this.f(r, prev[i]))
						+ (halfG * (this.f(r, prev.wrapAt(i+1)) + this.f(r, prev.wrapAt(i-1))));
		});
		^next;
	}

	embedInStream { |inval|
		var items = list.copy;
		var rStrm = r.asStream, gStrm = g.asStream;
		var rVal, gVal;
		repeats.do({
			rVal = rStrm.next(inval);
			gVal = gStrm.next(inval);
			if (rVal.isNil || gVal.isNil, { ^inval });
			items = this.evolve(items, rVal, gVal);
			inval = items.yield;
		});
		^inval;
	}

	plot { |n=500,skip=0|
		var rct, cell, stream = this.asStream, strVal;
		var win;
		rct = Rect(0, 0, list.size, n);
		win = Window("r: " ++ r ++ " g: " ++ g ++ " size: " ++ list.size ++ "x" ++ n, rct, false);
		win.view.background = Color.black;
		skip.do { stream.next };
		win.drawFunc = {
			n.do {|i|
				strVal = stream.next;
				strVal.do {|item,j|
					Pen.fillColor = Color.gray(item);
					Pen.fillRect(Rect(j, i, 1, 1));
				};
			};

		};
		win.front;
		CmdPeriod.doOnce { win.close };
	}
}

Pgcm : Pcml {

	evolve { |prev, r, g|
		var next = Array.newClear(prev.size), nG = g / prev.size, sum = 0;
		prev.do { |item, i| sum = sum + this.f(r, item) };
		prev.do { |item, i| next[i] = ((1.0 - g) * this.f(r, item)) + (nG * sum) };
		^next;
	}
}
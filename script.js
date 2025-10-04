// Audio STFT demo: drag/drop WAV, 512-bin STFT -> iSTFT, playback, and waveform with playhead

(() => {
	// ---------- DOM ----------
	const dropZone = document.getElementById('dropZone');
	const fileInput = document.getElementById('fileInput');
	const browseBtn = document.getElementById('browseBtn');
	const statusEl = document.getElementById('status');
	const fileInfoEl = document.getElementById('fileInfo');
	const fileNameEl = document.getElementById('fileName');
	const fileDurationEl = document.getElementById('fileDuration');
	const fileSampleRateEl = document.getElementById('fileSampleRate');
	const controlsEl = document.getElementById('controls');
	const playBtn = document.getElementById('playBtn');
	const pauseBtn = document.getElementById('pauseBtn');
	const downloadBtn = document.getElementById('downloadBtn');
	const waveformContainerBefore = document.getElementById('waveformContainerBefore');
	const waveformCanvasBefore = document.getElementById('waveformCanvasBefore');
	const waveformContainer = document.getElementById('waveformContainer');
	const waveformCanvas = document.getElementById('waveformCanvas');
	// Spectrogram DOM
	const spectrogramContainerBefore = document.getElementById('spectrogramContainerBefore');
	const spectrogramCanvasBefore = document.getElementById('spectrogramCanvasBefore');
	const spectrogramContainerAfter = document.getElementById('spectrogramContainerAfter');
	const spectrogramCanvasAfter = document.getElementById('spectrogramCanvasAfter');

	// ---------- Audio state ----------
	let audioCtx = null;
	let srcNode = null; // Current BufferSource for playback
	let originalBuffer = null; // Decoded source
	let reconBuffer = null; // Reconstructed buffer after iSTFT
	let isPlaying = false;
	let playStartTime = 0; // audioCtx.currentTime when playback started
	let pauseOffset = 0; // seconds into the buffer when paused/seeked
		// Expose last STFT frames for future experiments
		let lastSTFT = null; // Array per-channel: [ [{re,im}, ...], ... ]

	// ---------- STFT config ----------
	const FFT_SIZE = 512; // bins
	const HOP_SIZE = FFT_SIZE / 2; // 50% overlap
	const hannWindow = createHannWindow(FFT_SIZE);

	// Canvas dimensions
	const CANVAS_HEIGHT = 150; // px (75% of previous 200)

	// Precompute FFT tables for performance
	const fftTables = createFFT(FFT_SIZE);

	// ---------- UI wiring ----------
	// Drag & drop
	;['dragenter', 'dragover'].forEach(evt => {
		dropZone.addEventListener(evt, e => {
			e.preventDefault();
			e.stopPropagation();
			dropZone.classList.add('drag-over');
		});
	});
	;['dragleave', 'drop'].forEach(evt => {
		dropZone.addEventListener(evt, e => {
			e.preventDefault();
			e.stopPropagation();
			if (evt === 'drop') {
				dropZone.classList.remove('drag-over');
				const files = e.dataTransfer?.files;
				if (files && files.length) handleFile(files[0]);
			} else {
				dropZone.classList.remove('drag-over');
			}
		});
	});
	// clicking the drop zone opens the file chooser; ensure the inner Browse button
	// doesn't bubble its click up to the drop zone (which would call fileInput.click() twice)
	dropZone.addEventListener('click', () => fileInput.click());
	if (browseBtn) {
		browseBtn.addEventListener('click', (e) => {
			e.stopPropagation();
			e.preventDefault();
			fileInput.click();
		});
	}
	fileInput.addEventListener('change', () => {
		const file = fileInput.files?.[0];
		if (file) handleFile(file);
	});

	// Playback controls
	playBtn.addEventListener('click', () => {
		if (reconBuffer) {
			play();
		}
	});
	pauseBtn.addEventListener('click', () => pause());
	downloadBtn?.addEventListener('click', async () => {
		if (!reconBuffer) return;
		try {
			const wavBlob = audioBufferToWavBlob(reconBuffer);
			const baseName = (fileNameEl.textContent || 'audio').replace(/\.[^.]+$/, '');
			const outName = `${baseName}-lindverted.wav`;
			triggerDownload(wavBlob, outName);
		} catch (e) {
			console.error(e);
			setStatus('Failed to export WAV', 'error');
		}
	});

	// Waveform interactions
	window.addEventListener('resize', debounce(() => {
		redrawWaveforms();
		redrawSpectrograms();
	}, 150));
	waveformCanvas.addEventListener('click', (e) => {
		if (!reconBuffer || !audioCtx) return;
		const rect = waveformCanvas.getBoundingClientRect();
		const x = e.clientX - rect.left;
		const ratio = clamp(x / rect.width, 0, 1);
		const newTime = ratio * reconBuffer.duration;
		const wasPlaying = isPlaying;
		stopPlaybackInternal();
		pauseOffset = newTime;
		if (wasPlaying) play(); else redrawWaveforms();
	});

	// Animation for playhead
	let rafId = 0;
	function startRAF() {
		cancelAnimationFrame(rafId);
		const loop = () => {
			redrawWaveforms();
			rafId = requestAnimationFrame(loop);
		};
		rafId = requestAnimationFrame(loop);
	}
	function stopRAF() {
		cancelAnimationFrame(rafId);
	}

	// ---------- Core flow ----------
	async function handleFile(file) {
		resetUI();
		setStatus(`Reading ${file.name} …`, 'processing');

		if (!audioCtx) audioCtx = new (window.AudioContext || window.webkitAudioContext)();

		try {
			const arrayBuffer = await readFileAsArrayBuffer(file);
			setStatus('Decoding audio … This may take a minute.', 'processing');
			const decoded = await audioCtx.decodeAudioData(arrayBuffer.slice(0));
			originalBuffer = decoded;

			// Show file info
			fileNameEl.textContent = file.name;
			fileDurationEl.textContent = formatTime(decoded.duration);
			fileSampleRateEl.textContent = `${decoded.sampleRate.toLocaleString()} Hz`;
			fileInfoEl.style.display = 'block';

			// Draw 'before' waveform immediately; 'after' once processed
			waveformContainerBefore.style.display = 'block';
			drawWaveformOnCanvas(waveformCanvasBefore, decoded.getChannelData(0), decoded.sampleRate, false);
			// Draw 'before' spectrogram immediately
			spectrogramContainerBefore.style.display = 'block';
			drawSpectrogramOnCanvas(spectrogramCanvasBefore, decoded.getChannelData(0), decoded.sampleRate);

			// STFT -> mag/phase (apply crude highpass) -> complex -> iSTFT
			// setStatus('Processing STFT → mag/phase (highpass) → complex → iSTFT (512 bins, 50% overlap) …', 'processing');
			const recon = await stftProcessBuffer(decoded);
			reconBuffer = recon;

			// Compare errors
			const { rms, peak } = compareBuffers(decoded, recon);
			const msg = `Done. Diff — RMS: ${rms.toExponential(2)}, Peak: ${peak.toExponential(2)}`;
			setStatus('');

			// Show controls and 'after' waveform (reconstructed)
			controlsEl.style.display = 'flex';
			downloadBtn.style.display = 'inline-flex';
			waveformContainer.style.display = 'block';
			redrawWaveforms();
			// Draw 'after' spectrogram
			spectrogramContainerAfter.style.display = 'block';
			drawSpectrogramOnCanvas(spectrogramCanvasAfter, recon.getChannelData(0), recon.sampleRate);
			startRAF();
		} catch (err) {
			console.error(err);
			setStatus(`Error: ${err?.message || err}`, 'error');
		}
	}

	function modifyInPlace(magPhaseFrames) {
		for (let f = 0; f < magPhaseFrames.length; f++) {
			const mag = magPhaseFrames[f].mag;
			const binCount = mag.length;

			// Crude high-pass example: zero magnitudes for bins with |k| <= cutoffBins, preserving conjugate symmetry.
			// cutoffFrac is 0..1 relative to the positive-frequency band (0..N/2).
			// COMMENTED OUT: do not delete
			// const cutoffFrac = 0.1; // fraction of positive frequencies to zero out
			// const N = mag.length;
			// const half = N >> 1; // N/2 (Nyquist index)
			// const cutoff = Math.max(0, Math.min(half, Math.floor(half * cutoffFrac)));
			// for (let k = 0; k <= cutoff; k++) mag[k] = 0;
			// for (let k = N - cutoff; k < N; k++) mag[k] = 0;

			// --- Normalized "inversion" ---
			// Step 1: gather statistics for this frame.
			// `frameMax` anchors our inversion around the loudest bin in the frame, and
			// `originalSum` lets us keep the frame's total magnitude (energy) unchanged later.
			let frameMax = 0;
			let originalSum = 0;
			for (let k = 0; k < binCount; k++) {
				const m = mag[k];
				if (m > frameMax) frameMax = m;
				originalSum += m;
			}
			// If the frame is effectively silent, leave it untouched to avoid noise injection.
			if (frameMax <= 1e-9) continue;
			// Step 2: invert each bin relative to the loudest bin.
			// We normalize magnitudes so values near the max map close to 1, then flip (1 - norm).
			// Multiplying by `frameMax` restores the original units so the spectral shape is mirrored
			// without forcing everything into the [0,1] range.
			const denom = frameMax + 1e-6; // stabilize division when magnitudes are tiny
			const inverted = new Float32Array(binCount);
			let invertedSum = 0;
			for (let k = 0; k < binCount; k++) {
				const norm = mag[k] / denom; // 0..~1 for bins near the frame max
				const invMag = (1 - norm) * frameMax;
				inverted[k] = invMag;
				invertedSum += invMag;
			}
			// Step 3: rescale the inverted frame so total energy matches the original sum.
			// This prevents wild gain jumps when many bins get boosted by the inversion.
			const scale = invertedSum > 0 ? originalSum / invertedSum : 0;
			for (let k = 0; k < binCount; k++) mag[k] = inverted[k] * scale;
			// Phase left intact; array reference retained so downstream conversions stay valid.
		}
	}

	async function stftProcessBuffer(buffer) {
		const { numberOfChannels, length, sampleRate } = buffer;
			const outChannels = [];
			const stftPerChannel = [];
			const magPhasePerChannel = [];
		for (let ch = 0; ch < numberOfChannels; ch++) {
			const input = buffer.getChannelData(ch);
			// Pad at head and tail to reduce boundary transients; later we trim back
			const pad = FFT_SIZE;
			const padded = new Float32Array(length + 2 * pad);
			padded.set(input, pad);
			const frames = stft(padded, FFT_SIZE, HOP_SIZE, hannWindow);
			stftPerChannel.push(frames);
			const magPhaseFrames = framesToMagPhase(frames);
			modifyInPlace(magPhaseFrames);
			magPhasePerChannel.push(magPhaseFrames);
			const complexFrames = magPhaseToFrames(magPhaseFrames);
			const fullOut = istft(complexFrames, FFT_SIZE, HOP_SIZE, hannWindow);
			// Trim padding: keep original range [pad, pad+length)
			const output = fullOut.subarray(pad, pad + length);
			outChannels.push(output);
		}
			lastSTFT = stftPerChannel;
			// Expose a minimal debug handle on window for experimentation later
			try {
				window.AudioSTFTDebug = {
					frames: lastSTFT,
					magphase: magPhasePerChannel,
					fftSize: FFT_SIZE,
					hopSize: HOP_SIZE,
					window: 'hann'
				};
			} catch {}
		// Create an AudioBuffer from reconstructed channels
		const outLen = outChannels[0].length;
		const outBuffer = audioCtx.createBuffer(numberOfChannels, outLen, sampleRate);
		for (let ch = 0; ch < numberOfChannels; ch++) {
			outBuffer.getChannelData(ch).set(outChannels[ch]);
		}
		return outBuffer;
	}

	// ---------- Playback ----------
	function play() {
		if (!audioCtx || !reconBuffer) return;
		if (isPlaying) return;
		stopPlaybackInternal();
		srcNode = audioCtx.createBufferSource();
		srcNode.buffer = reconBuffer;
		srcNode.connect(audioCtx.destination);
		const offset = clamp(pauseOffset, 0, reconBuffer.duration);
		const duration = reconBuffer.duration - offset;
		srcNode.start(0, offset, duration);
		playStartTime = audioCtx.currentTime - offset;
		isPlaying = true;
		playBtn.style.display = 'none';
		pauseBtn.style.display = 'inline-flex';
		srcNode.onended = () => {
			if (!isPlaying) return; // if we paused via stop(), ignore
			stopPlaybackInternal();
			pauseOffset = 0;
			redrawWaveforms();
		};
	}

	function pause() {
		if (!isPlaying) return;
		if (!audioCtx) return;
		pauseOffset = getCurrentTime();
		stopPlaybackInternal();
		redrawWaveforms();
	}

	function stopPlaybackInternal() {
		if (srcNode) {
			try { srcNode.stop(); } catch {}
			try { srcNode.disconnect(); } catch {}
		}
		srcNode = null;
		isPlaying = false;
		playBtn.style.display = 'inline-flex';
		pauseBtn.style.display = 'none';
	}

	function getCurrentTime() {
		if (!audioCtx || !reconBuffer) return 0;
		if (!isPlaying) return pauseOffset;
		const t = audioCtx.currentTime - playStartTime;
		return clamp(t, 0, reconBuffer.duration);
	}

	// ---------- STFT / iSTFT ----------
	function stft(signal, nfft, hop, window) {
		const frames = [];
		const win = window;
		const re = new Float32Array(nfft);
		const im = new Float32Array(nfft);
		for (let start = 0; start + nfft <= signal.length; start += hop) {
			// Copy windowed frame
			for (let i = 0; i < nfft; i++) {
				re[i] = signal[start + i] * win[i];
				im[i] = 0;
			}
			fftInPlace(re, im, fftTables);
			// Save a copy of this frame's spectrum
			frames.push({ re: new Float32Array(re), im: new Float32Array(im) });
		}
		// If tail remains, zero-pad final frame
		const remainder = signal.length % hop;
		if (signal.length >= nfft && remainder !== 0 && (signal.length - nfft) % hop !== 0) {
			const start = Math.max(0, signal.length - nfft);
			for (let i = 0; i < nfft; i++) {
				const s = start + i < signal.length ? signal[start + i] : 0;
				re[i] = s * win[i];
				im[i] = 0;
			}
			fftInPlace(re, im, fftTables);
			frames.push({ re: new Float32Array(re), im: new Float32Array(im) });
		}
		return frames;
	}

	function istft(frames, nfft, hop, window, targetLen) {
		if (!frames.length) return new Float32Array(0);
		const win = window;
		const outLen = (frames.length - 1) * hop + nfft;
		const out = new Float32Array(outLen);
		const norm = new Float32Array(outLen);

		const re = new Float32Array(nfft);
		const im = new Float32Array(nfft);
		for (let f = 0; f < frames.length; f++) {
			re.set(frames[f].re);
			im.set(frames[f].im);
			ifftInPlace(re, im, fftTables);
			// Overlap-add with synthesis window (same Hann)
			const start = f * hop;
			for (let i = 0; i < nfft; i++) {
				const w = win[i];
				out[start + i] += re[i] * w; // use time-domain frame times window
				norm[start + i] += w * w; // accumulate window-squared for COLA normalization
			}
		}

		// Normalize to correct for overlap
		for (let i = 0; i < outLen; i++) {
			const d = norm[i];
			if (d > 1e-8) out[i] /= d;
		}

		// Trim or pad to original length
		if (typeof targetLen === 'number') {
			if (outLen > targetLen) return out.subarray(0, targetLen);
			if (outLen < targetLen) {
				const padded = new Float32Array(targetLen);
				padded.set(out);
				return padded;
			}
		}
		return out;
	}

	// ---------- Complex <-> Mag/Phase helpers ----------
	function framesToMagPhase(frames) {
		const out = new Array(frames.length);
		for (let f = 0; f < frames.length; f++) {
			const re = frames[f].re;
			const im = frames[f].im;
			const N = re.length;
			const mag = new Float32Array(N);
			const pha = new Float32Array(N);
			for (let k = 0; k < N; k++) {
				const r = re[k];
				const ii = im[k];
				mag[k] = Math.hypot(r, ii);
				pha[k] = Math.atan2(ii, r);
			}
			out[f] = { mag, pha };
		}
		return out;
	}

	function magPhaseToFrames(frames) {
		const out = new Array(frames.length);
		for (let f = 0; f < frames.length; f++) {
			const mag = frames[f].mag;
			const pha = frames[f].pha;
			const N = mag.length;
			const re = new Float32Array(N);
			const im = new Float32Array(N);
			for (let k = 0; k < N; k++) {
				const m = mag[k];
				const p = pha[k];
				re[k] = m * Math.cos(p);
				im[k] = m * Math.sin(p);
			}
			out[f] = { re, im };
		}
		return out;
	}

	// ---------- FFT implementation (radix-2, in-place) ----------
	function createFFT(N) {
		if ((N & (N - 1)) !== 0) throw new Error('FFT size must be power of 2');
		const bitCount = Math.log2(N) | 0;
		const bitRev = new Uint32Array(N);
		for (let i = 0; i < N; i++) bitRev[i] = reverseBits(i, bitCount);
		const cos = new Float32Array(N / 2);
		const sin = new Float32Array(N / 2);
		for (let k = 0; k < N / 2; k++) {
			const ang = -2 * Math.PI * k / N;
			cos[k] = Math.cos(ang);
			sin[k] = Math.sin(ang);
		}
		return { N, bitCount, bitRev, cos, sin };
	}

	function fftInPlace(re, im, tbl) {
		const { N, bitRev, cos, sin } = tbl;
		// Bit-reversal permutation
		for (let i = 0; i < N; i++) {
			const j = bitRev[i];
			if (j > i) {
				const tr = re[i]; re[i] = re[j]; re[j] = tr;
				const ti = im[i]; im[i] = im[j]; im[j] = ti;
			}
		}
		// Cooley–Tukey
		for (let size = 2; size <= N; size <<= 1) {
			const half = size >> 1;
			const step = N / size;
			for (let i = 0; i < N; i += size) {
				for (let j = 0, k = 0; j < half; j++, k += step) {
					const tcos = cos[k];
					const tsin = sin[k];
					const ur = re[i + j];
					const ui = im[i + j];
					const vr = re[i + j + half];
					const vi = im[i + j + half];
					// t = W * v
					const tr = tcos * vr - tsin * vi;
					const ti = tcos * vi + tsin * vr;
					re[i + j] = ur + tr;
					im[i + j] = ui + ti;
					re[i + j + half] = ur - tr;
					im[i + j + half] = ui - ti;
				}
			}
		}
	}

	function ifftInPlace(re, im, tbl) {
		// Conjugate
		for (let i = 0; i < tbl.N; i++) im[i] = -im[i];
		fftInPlace(re, im, tbl);
		// Conjugate again and scale
		const invN = 1 / tbl.N;
		for (let i = 0; i < tbl.N; i++) {
			re[i] = re[i] * invN; // im[i] would be -im[i]*invN but we don't need im after iFFT
			im[i] = -im[i] * invN;
		}
	}

	function reverseBits(x, bits) {
		let y = 0;
		for (let i = 0; i < bits; i++) {
			y = (y << 1) | (x & 1);
			x >>= 1;
		}
		return y >>> 0;
	}

	function createHannWindow(N) {
		const w = new Float32Array(N);
		for (let n = 0; n < N; n++) w[n] = 0.5 * (1 - Math.cos(2 * Math.PI * n / (N - 1)));
		return w;
	}

	// ---------- Waveform drawing ----------
	function drawWaveform(signal, sampleRate) {
		const dpr = window.devicePixelRatio || 1;
		const cssWidth = waveformContainer.clientWidth || 600;
		const cssHeight = CANVAS_HEIGHT;
		waveformCanvas.width = Math.floor(cssWidth * dpr);
		waveformCanvas.height = Math.floor(cssHeight * dpr);
		waveformCanvas.style.width = cssWidth + 'px';
		waveformCanvas.style.height = cssHeight + 'px';
		const ctx = waveformCanvas.getContext('2d');
		if (!ctx) return;
		ctx.scale(dpr, dpr);
		ctx.clearRect(0, 0, cssWidth, cssHeight);

		// Background
		ctx.fillStyle = '#eaeaff';
		ctx.fillRect(0, 0, cssWidth, cssHeight);

		// Compute min/max per pixel column
		const pixels = cssWidth;
		const step = Math.max(1, Math.floor(signal.length / pixels));
		const halfH = cssHeight / 2;
		ctx.strokeStyle = '#667eea';
		ctx.lineWidth = 1;
		ctx.beginPath();
		for (let x = 0; x < pixels; x++) {
			const start = x * step;
			let min = 1, max = -1;
			for (let i = 0; i < step && start + i < signal.length; i++) {
				const v = signal[start + i];
				if (v < min) min = v;
				if (v > max) max = v;
			}
			const y1 = halfH - min * halfH;
			const y2 = halfH - max * halfH;
			ctx.moveTo(x + 0.5, y1);
			ctx.lineTo(x + 0.5, y2);
		}
		ctx.stroke();

		// Playhead only on 'after' canvas
		if (reconBuffer) {
			const t = getCurrentTime();
			const ratio = reconBuffer.duration ? t / reconBuffer.duration : 0;
			const x = Math.round(ratio * cssWidth) + 0.5;
			ctx.strokeStyle = '#f59e0b';
			ctx.lineWidth = 2;
			ctx.beginPath();
			ctx.moveTo(x, 0);
			ctx.lineTo(x, cssHeight);
			ctx.stroke();
		}
	}

	function drawWaveformOnCanvas(canvas, signal, sampleRate, drawPlayhead) {
		const dpr = window.devicePixelRatio || 1;
		const container = canvas.parentElement;
		const cssWidth = (container?.clientWidth) || 600;
		const cssHeight = CANVAS_HEIGHT;
		canvas.width = Math.floor(cssWidth * dpr);
		canvas.height = Math.floor(cssHeight * dpr);
		canvas.style.width = cssWidth + 'px';
		canvas.style.height = cssHeight + 'px';
		const ctx = canvas.getContext('2d');
		if (!ctx) return;
		ctx.setTransform(1, 0, 0, 1, 0, 0);
		ctx.scale(dpr, dpr);
		ctx.clearRect(0, 0, cssWidth, cssHeight);

		ctx.fillStyle = '#eaeaff';
		ctx.fillRect(0, 0, cssWidth, cssHeight);

		const pixels = cssWidth;
		const step = Math.max(1, Math.floor(signal.length / pixels));
		const halfH = cssHeight / 2;
		ctx.strokeStyle = '#667eea';
		ctx.lineWidth = 1;
		ctx.beginPath();
		for (let x = 0; x < pixels; x++) {
			const start = x * step;
			let min = 1, max = -1;
			for (let i = 0; i < step && start + i < signal.length; i++) {
				const v = signal[start + i];
				if (v < min) min = v;
				if (v > max) max = v;
			}
			const y1 = halfH - min * halfH;
			const y2 = halfH - max * halfH;
			ctx.moveTo(x + 0.5, y1);
			ctx.lineTo(x + 0.5, y2);
		}
		ctx.stroke();

		if (drawPlayhead && reconBuffer) {
			const t = getCurrentTime();
			const ratio = reconBuffer.duration ? t / reconBuffer.duration : 0;
			const x = Math.round(ratio * cssWidth) + 0.5;
			ctx.strokeStyle = '#f59e0b';
			ctx.lineWidth = 2;
			ctx.beginPath();
			ctx.moveTo(x, 0);
			ctx.lineTo(x, cssHeight);
			ctx.stroke();
		}
	}

	function redrawWaveforms() {
		if (originalBuffer) {
			waveformContainerBefore.style.display = 'block';
			drawWaveformOnCanvas(waveformCanvasBefore, originalBuffer.getChannelData(0), originalBuffer.sampleRate, false);
		}
		if (reconBuffer) {
			waveformContainer.style.display = 'block';
			drawWaveformOnCanvas(waveformCanvas, reconBuffer.getChannelData(0), reconBuffer.sampleRate, true);
		}
	}

	// ---------- Spectrogram drawing ----------
	function drawSpectrogramOnCanvas(canvas, signal, sampleRate) {
		// Use the same STFT parameters for the image to keep it consistent
		const frames = stft(signal, FFT_SIZE, HOP_SIZE, hannWindow);
			// Use only the first half (including Nyquist) for real signals
			const bins = frames.length ? frames[0].re.length : FFT_SIZE;
			const half = (bins >> 1); // N/2 index (Nyquist)
			const usedBins = half + 1; // 0..N/2 inclusive
			const mags = frames.map(fr => {
				const m = new Float32Array(usedBins);
				for (let k = 0; k < usedBins; k++) m[k] = Math.hypot(fr.re[k], fr.im[k]);
				return m;
			});
			// Convert to log scale and normalize per full image (only used bins)
			let max = 1e-12;
			for (const m of mags) for (let k = 0; k < m.length; k++) if (m[k] > max) max = m[k];
			const logMags = mags.map(m => {
				const out = new Float32Array(m.length);
				for (let k = 0; k < m.length; k++) out[k] = Math.log10(1e-12 + m[k] / max);
				return out;
			});

		const dpr = window.devicePixelRatio || 1;
		const container = canvas.parentElement;
		const cssWidth = (container?.clientWidth) || 600;
		const cssHeight = CANVAS_HEIGHT;
		canvas.width = Math.floor(cssWidth * dpr);
		canvas.height = Math.floor(cssHeight * dpr);
		canvas.style.width = cssWidth + 'px';
		canvas.style.height = cssHeight + 'px';
		const ctx = canvas.getContext('2d');
		if (!ctx) return;
		// Work directly in device pixels; do not scale after setting size
		ctx.setTransform(1,0,0,1,0,0);
		ctx.clearRect(0,0,canvas.width,canvas.height);
		// White background (match waveform)
		ctx.fillStyle = '#ffffff';
		ctx.fillRect(0,0,canvas.width,canvas.height);

		const cols = logMags.length;
		const rows = usedBins; // show 0..N/2
		const img = ctx.createImageData(canvas.width, canvas.height);
		// Purple colormap on white background
		for (let x = 0; x < canvas.width; x++) {
			const tIdx = Math.floor(x * cols / canvas.width);
			const col = logMags[Math.min(tIdx, cols - 1)];
			// Log-frequency mapping along vertical axis
			const fNyquist = sampleRate / 2;
			const fPerBin = fNyquist / half;
			const fMin = Math.max(20, fPerBin); // avoid log(0), clamp to at least one bin
			const fMax = fNyquist;
			const logBase = Math.log(fMax / fMin + 1e-12);
			for (let y = 0; y < canvas.height; y++) {
				const p = 1 - (y / (canvas.height - 1)); // top->1 (high), bottom->0 (low)
				const f = fMin * Math.exp(p * logBase); // logarithmic spacing
				let freqIdx = Math.round(f / fPerBin);
				if (!isFinite(freqIdx)) freqIdx = 0;
				freqIdx = Math.max(0, Math.min(rows - 1, freqIdx));
				const v = col[Math.min(freqIdx, col.length - 1)];
				const val = Math.max(0, Math.min(1, (v + 6) / 6)); // 0..1
				const purple = lerpColor([232, 234, 255], [102, 126, 234], val); // light -> brand purple
				const idx = (y * canvas.width + x) * 4;
				img.data[idx + 0] = purple[0];
				img.data[idx + 1] = purple[1];
				img.data[idx + 2] = purple[2];
				img.data[idx + 3] = 255;
			}
		}
		ctx.putImageData(img, 0, 0);
	}

	function lerpColor(c0, c1, t) {
		const r = Math.round(c0[0] + (c1[0] - c0[0]) * t);
		const g = Math.round(c0[1] + (c1[1] - c0[1]) * t);
		const b = Math.round(c0[2] + (c1[2] - c0[2]) * t);
		return [r, g, b];
	}

	function redrawSpectrograms() {
		if (originalBuffer) {
			spectrogramContainerBefore.style.display = 'block';
			drawSpectrogramOnCanvas(spectrogramCanvasBefore, originalBuffer.getChannelData(0), originalBuffer.sampleRate);
		}
		if (reconBuffer) {
			spectrogramContainerAfter.style.display = 'block';
			drawSpectrogramOnCanvas(spectrogramCanvasAfter, reconBuffer.getChannelData(0), reconBuffer.sampleRate);
		}
	}

	// ---------- Utils ----------
	function readFileAsArrayBuffer(file) {
		return new Promise((resolve, reject) => {
			const fr = new FileReader();
			fr.onload = () => resolve(fr.result);
			fr.onerror = () => reject(fr.error || new Error('Failed to read file'));
			fr.readAsArrayBuffer(file);
		});
	}

	function setStatus(text, type = '') {
		statusEl.textContent = text;
		statusEl.classList.remove('processing', 'error', 'success');
		if (type) statusEl.classList.add(type);
	}

	function resetUI() {
		stopRAF();
		stopPlaybackInternal();
		pauseOffset = 0;
		reconBuffer = null;
		originalBuffer = null;
		fileInfoEl.style.display = 'none';
		controlsEl.style.display = 'none';
		waveformContainer.style.display = 'none';
		setStatus('');
	}

	function compareBuffers(a, b) {
		const len = Math.min(a.length, b.length);
		const chs = Math.min(a.numberOfChannels, b.numberOfChannels);
		let sumSq = 0;
		let peak = 0;
		for (let ch = 0; ch < chs; ch++) {
			const A = a.getChannelData(ch);
			const B = b.getChannelData(ch);
			for (let i = 0; i < len; i++) {
				const d = A[i] - B[i];
				sumSq += d * d;
				const ad = Math.abs(d);
				if (ad > peak) peak = ad;
			}
		}
		const rms = Math.sqrt(sumSq / (len * chs));
		return { rms, peak };
	}

	function formatTime(sec) {
		const s = Math.floor(sec % 60).toString().padStart(2, '0');
		const m = Math.floor(sec / 60).toString();
		return `${m}:${s}`;
	}

	function clamp(x, a, b) { return Math.min(b, Math.max(a, x)); }
	function debounce(fn, ms) {
		let t; return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
	}

	// ---------- WAV export ----------
	function audioBufferToWavBlob(buffer) {
		const numChannels = buffer.numberOfChannels;
		const sampleRate = buffer.sampleRate;
		const length = buffer.length;
		// 16-bit PCM WAV
		const bytesPerSample = 2;
		const blockAlign = numChannels * bytesPerSample;
		const byteRate = sampleRate * blockAlign;
		const dataSize = length * blockAlign;
		const headerSize = 44;
		const totalSize = headerSize + dataSize;
		const ab = new ArrayBuffer(totalSize);
		const view = new DataView(ab);

		// Write WAV header (RIFF)
		writeString(view, 0, 'RIFF');
		view.setUint32(4, totalSize - 8, true);
		writeString(view, 8, 'WAVE');
		writeString(view, 12, 'fmt ');
		view.setUint32(16, 16, true); // PCM chunk size
		view.setUint16(20, 1, true);  // PCM format
		view.setUint16(22, numChannels, true);
		view.setUint32(24, sampleRate, true);
		view.setUint32(28, byteRate, true);
		view.setUint16(32, blockAlign, true);
		view.setUint16(34, bytesPerSample * 8, true);
		writeString(view, 36, 'data');
		view.setUint32(40, dataSize, true);

		// Interleave channels and write PCM
		const channels = [];
		for (let ch = 0; ch < numChannels; ch++) channels.push(buffer.getChannelData(ch));
		let offset = headerSize;
		for (let i = 0; i < length; i++) {
			for (let ch = 0; ch < numChannels; ch++) {
				const sample = Math.max(-1, Math.min(1, channels[ch][i]));
				view.setInt16(offset, sample < 0 ? sample * 0x8000 : sample * 0x7FFF, true);
				offset += 2;
			}
		}

		return new Blob([ab], { type: 'audio/wav' });
	}

	function writeString(view, offset, str) {
		for (let i = 0; i < str.length; i++) view.setUint8(offset + i, str.charCodeAt(i));
	}

	function triggerDownload(blob, filename) {
		const url = URL.createObjectURL(blob);
		const a = document.createElement('a');
		a.href = url;
		a.download = filename;
		document.body.appendChild(a);
		a.click();
		document.body.removeChild(a);
		URL.revokeObjectURL(url);
	}
})();


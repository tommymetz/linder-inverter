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
	const waveformContainer = document.getElementById('waveformContainer');
	const waveformCanvas = document.getElementById('waveformCanvas');

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
	dropZone.addEventListener('click', () => fileInput.click());
	browseBtn?.addEventListener('click', () => fileInput.click());
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

	// Waveform interactions
	window.addEventListener('resize', debounce(() => {
		if (reconBuffer) drawWaveform(reconBuffer.getChannelData(0), reconBuffer.sampleRate);
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
		if (wasPlaying) play(); else drawWaveform(reconBuffer.getChannelData(0), reconBuffer.sampleRate);
	});

	// Animation for playhead
	let rafId = 0;
	function startRAF() {
		cancelAnimationFrame(rafId);
		const loop = () => {
			if (reconBuffer) drawWaveform(reconBuffer.getChannelData(0), reconBuffer.sampleRate);
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
			setStatus('Decoding audio …', 'processing');
			const decoded = await audioCtx.decodeAudioData(arrayBuffer.slice(0));
			originalBuffer = decoded;

			// Show file info
			fileNameEl.textContent = file.name;
			fileDurationEl.textContent = formatTime(decoded.duration);
			fileSampleRateEl.textContent = `${decoded.sampleRate.toLocaleString()} Hz`;
			fileInfoEl.style.display = 'block';

			// Draw original waveform now (we will overwrite with reconstructed after)
			waveformContainer.style.display = 'block';
			drawWaveform(decoded.getChannelData(0), decoded.sampleRate);

			// STFT -> mag/phase (apply crude highpass) -> complex -> iSTFT
			setStatus('Processing STFT → mag/phase (highpass) → complex → iSTFT (512 bins, 50% overlap) …', 'processing');
			const recon = await stftProcessBuffer(decoded);
			reconBuffer = recon;

			// Compare errors
			const { rms, peak } = compareBuffers(decoded, recon);
			const msg = `Done. Applied crude highpass in mag/phase. Difference from original is expected. Error — RMS: ${rms.toExponential(2)}, Peak: ${peak.toExponential(2)}`;
			setStatus(msg);

			// Show controls and final waveform (reconstructed)
			controlsEl.style.display = 'flex';
			drawWaveform(recon.getChannelData(0), recon.sampleRate);
			startRAF();
		} catch (err) {
			console.error(err);
			setStatus(`Error: ${err?.message || err}`, 'error');
		}
	}

	function modifyInPlace(magPhaseFrames) {
		for (let f = 0; f < magPhaseFrames.length; f++) {
			const mag = magPhaseFrames[f].mag;

      // Crude high-pass: zero magnitudes for bins with |k| <= cutoffBins, preserving conjugate symmetry
      // cutoffFrac is 0..1 relative to the positive-frequency band (0..N/2). Example: 0.25 removes the lowest 25% of that band.
      const cutoffFrac = 0.25
			const N = mag.length;
			const half = N >> 1; // N/2 (Nyquist index)
			const cutoff = Math.max(0, Math.min(half, Math.floor(half * cutoffFrac)));
			for (let k = 0; k <= cutoff; k++) mag[k] = 0;
			for (let k = N - cutoff; k < N; k++) mag[k] = 0;
		}
	}

	async function stftProcessBuffer(buffer) {
		const { numberOfChannels, length, sampleRate } = buffer;
			const outChannels = [];
			const stftPerChannel = [];
			const magPhasePerChannel = [];
		for (let ch = 0; ch < numberOfChannels; ch++) {
			const input = buffer.getChannelData(ch);
			const frames = stft(input, FFT_SIZE, HOP_SIZE, hannWindow);
			stftPerChannel.push(frames);
			const magPhaseFrames = framesToMagPhase(frames);
			modifyInPlace(magPhaseFrames);
			magPhasePerChannel.push(magPhaseFrames);
			const complexFrames = magPhaseToFrames(magPhaseFrames);
			const output = istft(complexFrames, FFT_SIZE, HOP_SIZE, hannWindow, length);
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
			drawWaveform(reconBuffer.getChannelData(0), reconBuffer.sampleRate);
		};
	}

	function pause() {
		if (!isPlaying) return;
		if (!audioCtx) return;
		pauseOffset = getCurrentTime();
		stopPlaybackInternal();
		drawWaveform(reconBuffer.getChannelData(0), reconBuffer.sampleRate);
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
		const cssHeight = 200;
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

		// Playhead
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
})();


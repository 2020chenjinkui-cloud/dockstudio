import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  ApiError,
  createJob,
  createProject,
  fetchFile,
  getJob,
  listProjects,
} from './api';
import type { DockRequest } from './types';

function jsonResponse(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { 'Content-Type': 'application/json' },
  });
}

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('api client', () => {
  it('createProject 发送 JSON body 并返回 Project', async () => {
    const project = { id: 'p1', name: 'demo', path: '/data/p1', created_at: '2026-10-06' };
    const spy = vi.fn().mockResolvedValue(jsonResponse(project));
    vi.stubGlobal('fetch', spy);

    const res = await createProject('demo');
    expect(res).toEqual(project);

    const [url, init] = spy.mock.calls[0] as [string, RequestInit];
    expect(url).toBe('/api/projects');
    expect(init.method).toBe('POST');
    expect(JSON.parse(init.body as string)).toEqual({ name: 'demo' });
  });

  it('listProjects 返回数组', async () => {
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue(jsonResponse([])));
    await expect(listProjects()).resolves.toEqual([]);
  });

  it('getJob 命中 {detail} 错误时抛出 detail 文案', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn().mockResolvedValue(jsonResponse({ detail: '任务不存在' }, 404)),
    );
    const err = await getJob('missing').catch((e: unknown) => e);
    expect(err).toBeInstanceOf(ApiError);
    expect((err as ApiError).message).toBe('任务不存在');
    expect((err as ApiError).status).toBe(404);
  });

  it('非 JSON 错误体回退到状态行', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn().mockResolvedValue(new Response('oops', { status: 500, statusText: 'Internal Server Error' })),
    );
    const err = await getJob('x').catch((e: unknown) => e);
    expect((err as ApiError).message).toContain('500');
  });

  it('网络异常映射为"后端未启动"', async () => {
    vi.stubGlobal('fetch', vi.fn().mockRejectedValue(new TypeError('fetch failed')));
    const err = await listProjects().catch((e: unknown) => e);
    expect(err).toBeInstanceOf(ApiError);
    expect((err as ApiError).message).toContain('后端未启动');
    expect((err as ApiError).status).toBe(0);
  });

  it('createJob POST 到项目路径并携带对接请求', async () => {
    const spy = vi.fn().mockResolvedValue(jsonResponse({ job_id: 'j1' }));
    vi.stubGlobal('fetch', spy);
    const body: DockRequest = {
      engine: 'vina',
      receptor_id: 'r1',
      ligand_id: 'l1',
      box: { center: [1, 2, 3], size: [22, 22, 22] },
      params: {
        exhaustiveness: 8,
        num_modes: 9,
        energy_range: 3,
        cpu: 1,
        seed: null,
        scoring: 'vina',
      },
    };
    const res = await createJob('p1', body);
    expect(res.job_id).toBe('j1');
    const [url, init] = spy.mock.calls[0] as [string, RequestInit];
    expect(url).toBe('/api/projects/p1/jobs');
    expect(JSON.parse(init.body as string)).toEqual(body);
  });

  it('fetchFile 对路径做 URL 编码并返回纯文本', async () => {
    const spy = vi
      .fn()
      .mockResolvedValue(new Response('ATOM ...', { status: 200 }));
    vi.stubGlobal('fetch', spy);
    const text = await fetchFile('C:\\data\\a b.pdb');
    expect(text).toBe('ATOM ...');
    const [url] = spy.mock.calls[0] as [string];
    expect(url).toBe(`/api/file?path=${encodeURIComponent('C:\\data\\a b.pdb')}`);
  });
});

import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  ApiError,
  createJob,
  createProject,
  fetchFile,
  getJob,
  listProjects,
} from './api';
import type { BatchCreateRequest, DockRequest } from './types';

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

describe('api v0.3', () => {
  it('predictPockets POST receptor_id 与 top_n', async () => {
    const payload = {
      receptor_id: 'r1',
      method: 'builtin_grid',
      pockets: [
        { id: 'p1', center: [1, 2, 3], size: [20, 20, 20], score: 12.5, n_points: 42 },
      ],
    };
    const spy = vi.fn().mockResolvedValue(jsonResponse(payload));
    vi.stubGlobal('fetch', spy);
    const { predictPockets } = await import('./api');
    const res = await predictPockets('r1', 5);
    expect(res.pockets).toHaveLength(1);
    const [url, init] = spy.mock.calls[0] as [string, RequestInit];
    expect(url).toBe('/api/pockets');
    expect(JSON.parse(init.body as string)).toEqual({ receptor_id: 'r1', top_n: 5 });
  });

  it('createFlexReceptor 400 时透出 detail（含可用链）', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn().mockResolvedValue(jsonResponse({ detail: '残基不存在，可用链：A / B' }, 400)),
    );
    const { createFlexReceptor } = await import('./api');
    const err = await createFlexReceptor('r1', ['Z:1']).catch((e: unknown) => e);
    expect(err).toBeInstanceOf(ApiError);
    expect((err as ApiError).message).toBe('残基不存在，可用链：A / B');
    expect((err as ApiError).status).toBe(400);
  });

  it('getBatchResults 拼接 sort/order/limit/offset 查询串', async () => {
    const spy = vi.fn().mockResolvedValue(jsonResponse({ total: 0, rows: [] }));
    vi.stubGlobal('fetch', spy);
    const { getBatchResults } = await import('./api');
    await getBatchResults('b1', { sort: 'affinity', order: 'asc', limit: 50, offset: 100 });
    const [url] = spy.mock.calls[0] as [string];
    expect(url).toBe('/api/batches/b1/results?sort=affinity&order=asc&limit=50&offset=100');
  });

  it('getBatchResults 无参数时不带查询串', async () => {
    const spy = vi.fn().mockResolvedValue(jsonResponse({ total: 0, rows: [] }));
    vi.stubGlobal('fetch', spy);
    const { getBatchResults } = await import('./api');
    await getBatchResults('b1');
    const [url] = spy.mock.calls[0] as [string];
    expect(url).toBe('/api/batches/b1/results');
  });

  it('createBatch POST 到项目路径', async () => {
    const spy = vi.fn().mockResolvedValue(jsonResponse({ batch_id: 'b1' }));
    vi.stubGlobal('fetch', spy);
    const { createBatch } = await import('./api');
    const body: BatchCreateRequest = {
      engine: 'vina',
      receptor_id: 'r1',
      ligand_ids: ['l1', 'l2'],
      box: { center: [0, 0, 0], size: [22, 22, 22] },
      params: {
        exhaustiveness: 8,
        num_modes: 9,
        energy_range: 3,
        cpu: 1,
        seed: null,
        scoring: 'vina',
      },
    };
    const res = await createBatch('p1', body);
    expect(res.batch_id).toBe('b1');
    const [url, init] = spy.mock.calls[0] as [string, RequestInit];
    expect(url).toBe('/api/projects/p1/batches');
    expect(JSON.parse(init.body as string)).toEqual(body);
  });
});

describe('api v0.4', () => {
  it('buildReportQuery 拼接 sections 与 title', async () => {
    const { buildReportQuery } = await import('./api');
    expect(buildReportQuery({ sections: ['scores', 'diagram'], title: '我的报告' })).toBe(
      `sections=${encodeURIComponent('scores,diagram')}&title=${encodeURIComponent('我的报告')}`,
    );
  });

  it('buildReportQuery 空 sections / 空 title 时不带参数', async () => {
    const { buildReportQuery } = await import('./api');
    expect(buildReportQuery({ sections: [], title: '  ' })).toBe('');
    expect(buildReportQuery({})).toBe('');
  });

  it('reportUrlV2 组合路径与查询串', async () => {
    const { reportUrlV2 } = await import('./api');
    expect(reportUrlV2('j1', { sections: ['scores'] })).toBe('/api/jobs/j1/report?sections=scores');
    expect(reportUrlV2('j1')).toBe('/api/jobs/j1/report');
  });

  it('diagramUrl 拼 pose_index', async () => {
    const { diagramUrl } = await import('./api');
    expect(diagramUrl('j1', 3)).toBe('/api/jobs/j1/diagram?pose_index=3');
  });

  it('compareJobs 用逗号拼 job_ids', async () => {
    const spy = vi.fn().mockResolvedValue(jsonResponse({ rows: [] }));
    vi.stubGlobal('fetch', spy);
    const { compareJobs } = await import('./api');
    await compareJobs(['a', 'b', 'c']);
    const [url] = spy.mock.calls[0] as [string];
    expect(url).toBe('/api/compare?job_ids=a,b,c');
  });
});
